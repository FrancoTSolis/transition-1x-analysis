#!/usr/bin/env python3
"""Task 3: optimize=True labels run past the 500-iteration cap (same config, same canonical init, same optimizer).

This is generate_compressed_targets.py for one config (default square_reg0.005), with a passive L-BFGS-B callback
that records the objective after every iteration and snapshots the parameters at fixed iterations.  The callback does
not change the trajectory: where the run reproduces the stored label's floating-point conditions, the snapshot at
iteration 500 is bit-identical to the stored 500-iteration label (checked against --ref-dirs: ref_dU / ref_dZ).

Floating-point conditions matter.  The L-BFGS trajectory is chaotic: the XLA CPU thread pool (= the CPU affinity of
the process) changes the summation order, and a single-threaded run leaves a multi-threaded one after a few
iterations.  Pin every process with taskset (one core = the single-threaded path); unpinned, a process also uses
2-4 cores, because the XLA_FLAGS below do not stop XLA's own thread pool.

Output (one .npz per molecule), same layout and keys as generate_compressed_targets.py, so the directory works as a
--labels-dir (policy_dump.py, task3_build_tasks.py):
    <out-dir>/init/<name>.npz
    <out-dir>/<config>/<name>.npz   Z, U_re, U_im, dkappa_re, dkappa_im, resid, znorm, fun, nit, nfev, success, time
        + message, gmax (max |grad| at the end), maxiter, gtol, ftol (-1 = scipy default), trace_fun (objective
          after each iteration), snap_iters, snap_Z, snap_U_re, snap_U_im, snap_fun, snap_resid, snap_gmax,
          ref_dir, ref_nit, ref_dU, ref_dZ, ref_dresid (snapshot at the stored label's iteration vs that label)

Options beyond generate_compressed_targets.py: --gtol / --ftol (L-BFGS-B tolerances; scipy defaults 1e-5 / 2.2e-9),
--warm-dir (start from an existing label instead of the canonical init; fresh L-BFGS memory), --x0-noise / --x0-seed
(relative Gaussian perturbation of the initial parameter vector, e.g. 1e-15, to sample the trajectory noise).

    taskset -c <core> .lucj_venv/bin/python3 -u pretrain/followups/task3_fit_labels.py --names-file X \
        --out-dir rhf_targets_compressed_conv --maxiter 5000 --shard i --n-shards N --ref-dirs <stored label dirs>
"""
import argparse
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault(
    "XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402
import scipy.optimize  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from gauge_study.compressed_canonical import (  # noqa: E402
    canonical_exact_init, diag_coulomb_indices_for, rel_residual, relative_kappa, z_mask)


def parse_config(s: str):
    conn, reg = s.rsplit("_reg", 1)
    return conn, float(reg)


def compress_from_init_traced(t2, Z0, U0, *, connectivity="all-to-all", regularization=0.0, znorm_ref=None,
                              maxiter=500, snap_iters=(), method="L-BFGS-B", x64=True, gtol=None, ftol=None,
                              x0_noise=0.0, x0_seed=0):
    """Verbatim copy of gauge_study.compressed_canonical.compress_from_init plus a passive callback."""
    import jax
    import jax.numpy as jnp
    from opt_einsum import contract
    from ffsim.linalg.util import df_tensors_from_params, df_tensors_from_params_jax, \
        df_tensors_to_params

    if x64:
        jax.config.update("jax_enable_x64", True)

    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    n_tensors = Z0.shape[0]
    dci = diag_coulomb_indices_for(connectivity, norb)
    if znorm_ref is None:
        znorm_ref = float(np.sum(np.abs(Z0) ** 2))
    t2_j = jnp.asarray(t2)

    def fun(x):
        Z, U = df_tensors_from_params_jax(x, n_tensors, norb, dci)
        rec = 1j * contract("mpq,map,mip,mbq,mjq->ijab", Z, U, U.conj(), U, U.conj(),
                            optimize="greedy")[:nocc, :nocc, nocc:, nocc:]
        loss = 0.5 * jnp.sum(jnp.abs(rec - t2_j) ** 2)
        if regularization:
            loss += regularization * jnp.abs(jnp.sum(jnp.abs(Z) ** 2) - znorm_ref)
        return loss

    vg = jax.jit(jax.value_and_grad(fun))

    def vg_np(x):
        v, g = vg(jnp.asarray(x))
        return float(v), np.asarray(g, dtype=float)

    trace, snaps = [], {}
    want = set(int(s) for s in snap_iters)

    def callback(intermediate_result):          # new-style scipy callback: OptimizeResult(x, fun); passive
        trace.append(float(intermediate_result.fun))
        if len(trace) in want:
            snaps[len(trace)] = np.array(intermediate_result.x, dtype=float, copy=True)

    x0 = df_tensors_to_params(Z0, U0, dci)
    if x0_noise:                                 # round-off-level perturbation: samples the trajectory noise
        x0 = x0 * (1.0 + x0_noise * np.random.default_rng(x0_seed).standard_normal(x0.shape))
    t0 = time.time()
    opts = dict(maxiter=maxiter)                 # scipy defaults (ftol 2.2e-9, gtol 1e-5) unless given
    if gtol is not None:
        opts["gtol"] = gtol
    if ftol is not None:
        opts["ftol"] = ftol
    res = scipy.optimize.minimize(vg_np, x0, method=method, jac=True, options=opts, callback=callback)
    elapsed = time.time() - t0
    Z, U = df_tensors_from_params(res.x, n_tensors, norb, dci)

    def unpack(x):
        z, u = df_tensors_from_params(x, n_tensors, norb, dci)
        return np.asarray(z, dtype=float), np.asarray(u)

    snap_out = {k: (*unpack(x), trace[k - 1], float(np.max(np.abs(vg_np(x)[1])))) for k, x in sorted(snaps.items())}
    return dict(Z=np.asarray(Z, dtype=float), U=np.asarray(U), x=res.x,
                nit=int(res.nit), nfev=int(res.nfev), success=bool(res.success),
                message=str(res.message), fun=float(res.fun), time=elapsed,
                gmax=float(np.max(np.abs(vg_np(res.x)[1]))), trace=np.array(trace), snaps=snap_out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--names-file", nargs="+", required=True, help="one or more name files (concatenated)")
    ap.add_argument("--out-dir", default="rhf_targets_compressed_conv")
    ap.add_argument("--config", default="square_reg0.005")
    ap.add_argument("--maxiter", type=int, default=5000)
    ap.add_argument("--snap-iters", type=int, nargs="*", default=[500, 1000, 2000, 3000, 4000])
    ap.add_argument("--phase-mode", default="maxabs", choices=["maxabs", "diag"])
    ap.add_argument("--gtol", type=float, default=None, help="L-BFGS-B gtol (default: scipy's 1e-5)")
    ap.add_argument("--ftol", type=float, default=None, help="L-BFGS-B ftol (default: scipy's 2.2e-9)")
    ap.add_argument("--x0-noise", type=float, default=0.0,
                    help="relative Gaussian perturbation of the initial parameter vector (e.g. 1e-15)")
    ap.add_argument("--x0-seed", type=int, default=0)
    ap.add_argument("--warm-dir", default=None,
                    help="start L-BFGS from the label in <warm-dir>/<config>/<name>.npz instead of the canonical init "
                         "(fresh L-BFGS memory; the canonical init still sets the regularizer reference)")
    ap.add_argument("--ref-dirs", nargs="*", default=[],
                    help="existing label dirs; the first one holding <config>/<name>.npz is the reference")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    args = ap.parse_args()

    data_dir, out_dir = ROOT / args.data_dir, ROOT / args.out_dir
    (out_dir / "init").mkdir(parents=True, exist_ok=True)
    (out_dir / args.config).mkdir(parents=True, exist_ok=True)
    names = []
    for nf in args.names_file:
        names += [ln.strip() for ln in open(ROOT / nf) if ln.strip()]
    names = list(dict.fromkeys(names))[args.shard::args.n_shards]
    conn, reg = parse_config(args.config)

    t_start = time.time()
    n_done = n_skip = n_err = 0
    for i, name in enumerate(names):
        out_path = out_dir / args.config / f"{name}.npz"
        init_path = out_dir / "init" / f"{name}.npz"
        if out_path.exists() and init_path.exists():
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
            Z0, U0 = init.Z * z_mask(conn, norb)[None], init.U
            if args.warm_dir:
                w = np.load(ROOT / args.warm_dir / args.config / f"{name}.npz")
                Z0, U0 = w["Z"] * z_mask(conn, norb)[None], w["U_re"] + 1j * w["U_im"]
            r = compress_from_init_traced(
                t2, Z0, U0, connectivity=conn, regularization=reg,
                znorm_ref=init.znorm_full, maxiter=args.maxiter, snap_iters=args.snap_iters,
                gtol=args.gtol, ftol=args.ftol, x0_noise=args.x0_noise, x0_seed=args.x0_seed)
            dk = relative_kappa(init.U, r["U"])
            sn = r["snaps"]
            ks = sorted(sn)
            out = dict(
                Z=r["Z"], U_re=r["U"].real, U_im=r["U"].imag,
                dkappa_re=dk.real, dkappa_im=dk.imag,
                resid=rel_residual(t2, r["Z"], r["U"]),
                znorm=float(np.sum(r["Z"] ** 2)),
                fun=r["fun"], nit=r["nit"], nfev=r["nfev"],
                success=r["success"], time=r["time"], message=r["message"], gmax=r["gmax"],
                warm_dir=args.warm_dir or "", x0_noise=args.x0_noise, x0_seed=args.x0_seed, maxiter=args.maxiter, gtol=-1.0 if args.gtol is None else args.gtol,
                ftol=-1.0 if args.ftol is None else args.ftol, trace_fun=r["trace"],
                snap_iters=np.array(ks, dtype=int),
                snap_Z=np.array([sn[k][0] for k in ks]).reshape(len(ks), *r["Z"].shape),
                snap_U_re=np.array([sn[k][1].real for k in ks]).reshape(len(ks), *r["Z"].shape),
                snap_U_im=np.array([sn[k][1].imag for k in ks]).reshape(len(ks), *r["Z"].shape),
                snap_fun=np.array([sn[k][2] for k in ks]),
                snap_gmax=np.array([sn[k][3] for k in ks]),
                snap_resid=np.array([rel_residual(t2, sn[k][0], sn[k][1]) for k in ks]),
                nocc=nocc, nvirt=nvirt, norb=norb)
            # reproducibility check: snapshot 500 (or the end, if the fit stopped earlier) vs the existing label
            for rd in args.ref_dirs:
                rp = ROOT / rd / args.config / f"{name}.npz"
                if not rp.exists():
                    continue
                ref = np.load(rp)
                k_ref = int(ref["nit"])
                if k_ref in sn:
                    Zs, Us = sn[k_ref][0], sn[k_ref][1]
                elif k_ref == r["nit"]:
                    Zs, Us = r["Z"], r["U"]
                else:
                    break
                Ur = ref["U_re"] + 1j * ref["U_im"]
                out.update(ref_dir=rd, ref_nit=k_ref, ref_dU=float(np.max(np.abs(Us - Ur))),
                           ref_dZ=float(np.max(np.abs(Zs - ref["Z"]))),
                           ref_dresid=float(rel_residual(t2, Zs, Us) - float(ref["resid"])))
                break
            tmp = out_path.with_name(out_path.stem + ".tmp.npz")
            np.savez(tmp, **out)
            os.replace(tmp, out_path)
            n_done += 1
            print(f"  {name} norb {norb}: nit {r['nit']} nfev {r['nfev']} success {r['success']} "
                  f"resid {out['resid']:.5f} (snap {dict(zip(ks, np.round(out['snap_resid'], 5)))}) "
                  f"gmax {r['gmax']:.2e} {r['time']:.0f}s  [{r['message']}]"
                  + (f"  ref {out['ref_dir']}: dU {out['ref_dU']:.1e} dZ {out['ref_dZ']:.1e}" if "ref_dU" in out
                     else ""), flush=True)
        except Exception as e:  # noqa: BLE001
            n_err += 1
            print(f"  ERR {name}: {type(e).__name__}: {e}", flush=True)
        el = time.time() - t_start
        print(f"[shard {args.shard}] {i+1}/{len(names)} done={n_done} skip={n_skip} err={n_err}  {el:.0f}s",
              flush=True)
    print(f"[shard {args.shard}] DONE done={n_done} skip={n_skip} err={n_err} {time.time()-t_start:.0f}s",
          flush=True)


if __name__ == "__main__":
    main()
