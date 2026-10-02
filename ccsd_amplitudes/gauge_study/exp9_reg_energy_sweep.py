#!/usr/bin/env python3
"""Experiment 9: choose the compressed-DF regularization lambda by ENERGY.

Lin et al. (arXiv:2511.22476) add lambda*|sum||J||^2 - sum||J_exact||^2| to the
compressed-DF objective to stop the diagonal-Coulomb norm from blowing up
(Trotter error), and use lambda = 0.005 for N2.  One test molecule here showed
unregularized compression *destroys* the variational energy (-300% of the
correlation energy) while lambda=0.005 recovers 70% (vs 19.5% for the exact
n_reps=2 truncation).  This sweeps lambda on the smallest molecules of the
dataset (norb <= 16, exact statevector energies) and reports the fraction of
CCSD correlation energy recovered, the t2 residual and ||Z||.

Run from ccsd_amplitudes/ (ffsim venv):
    python3 -m gauge_study.exp9_reg_energy_sweep --n-procs 30
"""
from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402

from .compressed_canonical import canonical_exact_init, compress_from_init, rel_residual, z_mask  # noqa: E402

LAMBDAS = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]


_CORES: list[int] = []


def _pin_worker(cores):
    """Pin each pool worker to its own core: jax's CPU thread pool ignores
    intra_op_parallelism_threads, so unpinned workers oversubscribe the box."""
    import multiprocessing
    if not cores:
        return
    wid = multiprocessing.current_process()._identity[0] - 1
    try:
        os.sched_setaffinity(0, {cores[wid % len(cores)]})
    except Exception:  # noqa: BLE001
        pass


def job(args):
    name, data_dir, ham_dir, conn, lam, maxiter = args
    from pretrain.rl.energy import exact_energy, make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    d = np.load(Path(data_dir) / f"{name}.npz")
    t2 = d["t2"].astype(np.float64); t1 = d["t1"].astype(np.float64)
    nocc, _, nvirt, _ = t2.shape; norb = nocc + nvirt
    ham, norb_h, nelec, e_hf, e_ccsd = load_hamiltonian(ham_dir, name)
    assert norb_h == norb
    init = canonical_exact_init(t2)
    m = z_mask(conn, norb)[None]
    out = dict(name=name, lam=lam, conn=conn, e_hf=e_hf, e_ccsd=e_ccsd, norb=norb)
    t0 = time.time()
    if lam is None:   # exact-DF init energy only
        e = exact_energy(ham, norb, nelec, make_ucj_op(init.Z * m, init.U, conn, t1=t1))
        out.update(E=e, resid=rel_residual(t2, init.Z * m, init.U),
                   znorm=float(np.linalg.norm(init.Z * m)), nit=0)
    else:
        r = compress_from_init(t2, init.Z * m, init.U, connectivity=conn, regularization=lam,
                               znorm_ref=init.znorm_full, maxiter=maxiter)
        e = exact_energy(ham, norb, nelec, make_ucj_op(r["Z"], r["U"], conn, t1=t1))
        out.update(E=e, resid=rel_residual(t2, r["Z"], r["U"]), znorm=float(np.linalg.norm(r["Z"])),
                   nit=r["nit"])
    out["corr_frac"] = (e_hf - e) / (e_hf - e_ccsd)
    out["znorm_full"] = float(np.sqrt(init.znorm_full))
    out["t"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--ham-dir", default="rhf_hamiltonians")
    ap.add_argument("--names-file", default="gauge_study/names_small_norb_le16.txt")
    ap.add_argument("--connectivities", nargs="+", default=["square"])
    ap.add_argument("--lambdas", nargs="+", type=float, default=LAMBDAS)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--n-procs", type=int, default=30)
    ap.add_argument("--out", default="gauge_study/exp9_results.jsonl")
    args = ap.parse_args()
    names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    names = [n for n in names if (Path(args.ham_dir) / f"{n}.npz").exists()]
    jobs = [(n, args.data_dir, args.ham_dir, c, lam, args.maxiter)
            for n in names for c in args.connectivities for lam in [None] + list(args.lambdas)]
    # resume: skip (name, conn, lam) already in the jsonl
    res = []
    if Path(args.out).exists():
        for ln in open(args.out):
            try:
                res.append(json.loads(ln))
            except Exception:  # noqa: BLE001
                pass
    done = {(r["name"], r["conn"], r["lam"]) for r in res}
    jobs = [j for j in jobs if (j[0], j[3], j[4]) not in done]
    print(f"{len(names)} molecules x {len(args.connectivities)} conn x {len(args.lambdas)+1}: "
          f"{len(jobs)} jobs to run ({len(done)} already done)", flush=True)
    cores = sorted(os.sched_getaffinity(0))
    t0 = time.time()
    with Pool(args.n_procs, initializer=_pin_worker, initargs=(cores,)) as pool, open(args.out, "a") as f:
        for r in pool.imap_unordered(job, jobs):
            res.append(r); f.write(json.dumps(r) + "\n"); f.flush()
            print(f"  [{len(res)}/{len(jobs)}] {r['name']:22s} {r['conn']:10s} lam={r['lam']}  corr%={100*r['corr_frac']:7.1f}"
                  f"  resid={r['resid']:.3f} |Z|={r['znorm']:.2f} (full {r['znorm_full']:.2f}) nit={r['nit']} ({r['t']:.0f}s, {time.time()-t0:.0f}s)", flush=True)
    print(f"\n=== {len(names)} molecules (norb<=16), exact statevector energies; % of CCSD correlation energy recovered ===")
    print(f"{'conn':10s} {'lambda':>8s} {'median%':>8s} {'q25%':>7s} {'q75%':>7s} {'min%':>8s} {'resid':>6s} {'|Z|':>6s} {'nit':>5s}")
    for c in args.connectivities:
        for lam in [None] + list(args.lambdas):
            rs = [r for r in res if r["conn"] == c and r["lam"] == lam]
            if not rs:
                continue
            cf = 100 * np.array([r["corr_frac"] for r in rs])
            print(f"{c:10s} {('init' if lam is None else f'{lam:g}'):>8s} {np.median(cf):8.1f} {np.quantile(cf,.25):7.1f} "
                  f"{np.quantile(cf,.75):7.1f} {cf.min():8.1f} {np.median([r['resid'] for r in rs]):6.3f} "
                  f"{np.median([r['znorm'] for r in rs]):6.2f} {np.median([r['nit'] for r in rs]):5.0f}")
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
