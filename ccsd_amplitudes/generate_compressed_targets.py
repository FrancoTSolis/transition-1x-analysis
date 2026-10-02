#!/usr/bin/env python3
"""Generate canonical-gauge compressed (optimize=True) DF targets.

For each molecule in --data-dir (or the --names-file subset), computes the
canonical exact-DF init once and then, for every requested config
    <connectivity>_reg<regularization>
runs ffsim's compressed optimization (L-BFGS, same objective/regularizer as
ffsim.linalg.double_factorized_t2) from that canonical init.

Output layout (one .npz per molecule):
    <out-dir>/init/<name>.npz
        lam0, v0 (nocc,nvirt), gap, Z (2,n,n), U_re/U_im (2,n,n), w (2,n),
        znorm_full, resid_init_{all-to-all,square}
    <out-dir>/<conn>_reg<reg>/<name>.npz
        Z (2,n,n), U_re/U_im (2,n,n), dkappa_re/dkappa_im (2,n,n)
        [U = U_init expm(dkappa)], resid, znorm, nit, nfev, success, time
        [+ U2_re/U2_im, gauge_dist if --gauge-check: rerun from a random
           phase-rotated init; distance after phase re-alignment]

Shardable: --shard i --n-shards N.  Single-threaded BLAS/XLA per process.
"""
import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault(
    "XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gauge_study.compressed_canonical import (  # noqa: E402
    canonical_exact_init, compress_from_init, phase_dist, rel_residual,
    relative_kappa, z_mask)


def parse_config(s: str):
    conn, reg = s.rsplit("_reg", 1)
    return conn, float(reg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--names-file", default=None,
                    help="text file with one molecule name per line (subset)")
    ap.add_argument("--out-dir", default="rhf_targets_compressed")
    ap.add_argument("--configs", nargs="+",
                    default=["all-to-all_reg0", "all-to-all_reg0.005", "all-to-all_reg0.01",
                             "square_reg0", "square_reg0.005", "square_reg0.01"])
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--phase-mode", default="maxabs", choices=["maxabs", "diag"])
    ap.add_argument("--gauge-check", action="store_true",
                    help="also optimize from a random-phase-rotated init and "
                         "record the aligned distance (optimizer gauge sensitivity)")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--max-mols", type=int, default=0)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    (out_dir / "init").mkdir(parents=True, exist_ok=True)
    for c in args.configs:
        (out_dir / c).mkdir(parents=True, exist_ok=True)

    if args.names_file:
        names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    else:
        names = sorted(f.stem for f in data_dir.glob("*.npz") if f.stem != "_index")
    names = names[args.shard::args.n_shards]
    if args.max_mols:
        names = names[:args.max_mols]
    rng = np.random.default_rng(1234 + args.shard)

    t_start = time.time()
    n_done = n_skip = n_err = 0
    for i, name in enumerate(names):
        todo = [c for c in args.configs if not (out_dir / c / f"{name}.npz").exists()]
        init_path = out_dir / "init" / f"{name}.npz"
        if not todo and init_path.exists():
            n_skip += 1
            continue
        try:
            d = np.load(data_dir / f"{name}.npz")
            t2 = d["t2"].astype(np.float64)
            nocc, _, nvirt, _ = t2.shape
            norb = nocc + nvirt
            init = canonical_exact_init(t2, phase_mode=args.phase_mode)
            if not init_path.exists():
                np.savez(
                    init_path, lam0=init.lam0, v0=init.v0.reshape(nocc, nvirt),
                    gap=init.gap, Z=init.Z, U_re=init.U.real, U_im=init.U.imag,
                    w=init.w, znorm_full=init.znorm_full,
                    resid_init_all=rel_residual(t2, init.Z, init.U),
                    resid_init_square=rel_residual(
                        t2, init.Z * z_mask("square", norb)[None], init.U),
                    nocc=nocc, nvirt=nvirt, norb=norb)
            for c in todo:
                conn, reg = parse_config(c)
                Z0 = init.Z * z_mask(conn, norb)[None]
                r = compress_from_init(
                    t2, Z0, init.U, connectivity=conn, regularization=reg,
                    znorm_ref=init.znorm_full, maxiter=args.maxiter)
                dk = relative_kappa(init.U, r["U"])
                out = dict(
                    Z=r["Z"], U_re=r["U"].real, U_im=r["U"].imag,
                    dkappa_re=dk.real, dkappa_im=dk.imag,
                    resid=rel_residual(t2, r["Z"], r["U"]),
                    znorm=float(np.sum(r["Z"] ** 2)),
                    fun=r["fun"], nit=r["nit"], nfev=r["nfev"],
                    success=r["success"], time=r["time"],
                    nocc=nocc, nvirt=nvirt, norb=norb)
                if args.gauge_check:
                    ph = np.exp(2j * np.pi * rng.random((2, norb)))
                    U0b = init.U * ph[:, None, :]
                    r2 = compress_from_init(
                        t2, Z0, U0b, connectivity=conn, regularization=reg,
                        znorm_ref=init.znorm_full, maxiter=args.maxiter)
                    out.update(
                        U2_re=r2["U"].real, U2_im=r2["U"].imag, Z2=r2["Z"],
                        resid2=rel_residual(t2, r2["Z"], r2["U"]),
                        gauge_dist=float(np.sqrt(sum(
                            phase_dist(r["U"][k], r2["U"][k]) ** 2 for k in range(2)))),
                        gauge_dist_raw=float(np.linalg.norm(r["U"] - r2["U"])),
                        gauge_dZ=float(np.linalg.norm(r["Z"] - r2["Z"])))
                np.savez(out_dir / c / f"{name}.npz", **out)
            n_done += 1
        except Exception as e:  # noqa: BLE001
            n_err += 1
            print(f"  ERR {name}: {type(e).__name__}: {e}", flush=True)
        if (i + 1) % 10 == 0 or i == len(names) - 1:
            el = time.time() - t_start
            print(f"[shard {args.shard}] {i+1}/{len(names)} done={n_done} skip={n_skip} "
                  f"err={n_err}  {el:.0f}s  ({el/max(n_done,1):.1f}s/mol)", flush=True)
    print(f"[shard {args.shard}] DONE done={n_done} skip={n_skip} err={n_err} "
          f"{time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
