#!/usr/bin/env python3
"""Experiment 12: do the many compressed-DF minima differ in ENERGY?

exp10 showed every restart of the compressed optimization lands in a different
minimum (pairwise invariant-content distance ~0.2 ||t2||) at nearly the same
residual.  If those minima also share the same energy, the multimodality is
harmless (any of them is a lossless target); if not, energy -- not the residual
-- has to pick the point, which is the RL argument.

norb<=16 molecules only (exact statevector energies).  Per molecule: R restarts
(random phases + kappa kick) of the square, lambda=0.005 optimization, then the
exact LUCJ energy of each solution and of the canonical init.

Run from ccsd_amplitudes/ (ffsim venv):
    python3 -m gauge_study.exp12_valley_energy --names C2H3N_rxn2857_TS C2H3N_rxn2858_P --restarts 10 --n-procs 20
"""
from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402
import scipy.linalg  # noqa: E402

from .compressed_canonical import canonical_exact_init, compress_from_init, reconstruct_t2, rel_residual, z_mask  # noqa: E402


def _pin(cores):
    import multiprocessing
    if cores:
        wid = multiprocessing.current_process()._identity[0] - 1
        try:
            os.sched_setaffinity(0, {cores[wid % len(cores)]})
        except Exception:  # noqa: BLE001
            pass


def job(args):
    name, data_dir, ham_dir, seed, kick, lam, maxiter = args
    from pretrain.rl.energy import exact_energy, make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    d = np.load(Path(data_dir) / f"{name}.npz")
    t2 = d["t2"].astype(np.float64); t1 = d["t1"].astype(np.float64)
    nocc, _, nvirt, _ = t2.shape; norb = nocc + nvirt
    ham, _, nelec, e_hf, e_ccsd = load_hamiltonian(ham_dir, name)
    init = canonical_exact_init(t2)
    m = z_mask("square", norb)[None]
    rng = np.random.default_rng(max(seed, 0))
    t0 = time.time()
    if seed < 0:  # canonical init itself
        Z, U = init.Z * m, init.U
        nit = 0
    else:
        U0 = init.U.copy()
        if seed > 0:
            ph = np.exp(2j * np.pi * rng.random((2, norb)))
            U0 = U0 * ph[:, None, :]
            for k in range(2):
                A = rng.standard_normal((norb, norb)) + 1j * rng.standard_normal((norb, norb))
                U0[k] = U0[k] @ scipy.linalg.expm(kick * (A - A.conj().T) / np.sqrt(2 * norb))
        r = compress_from_init(t2, init.Z * m, U0, connectivity="square", regularization=lam,
                               znorm_ref=init.znorm_full, maxiter=maxiter)
        Z, U, nit = r["Z"], r["U"], r["nit"]
    e = exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1=t1))
    return dict(name=name, seed=seed, resid=rel_residual(t2, Z, U), znorm=float(np.linalg.norm(Z)),
                E=e, corr_frac=(e_hf - e) / (e_hf - e_ccsd), nit=nit, t=time.time() - t0,
                t2hat=reconstruct_t2(Z, U, nocc).ravel().tolist(), t2norm=float(np.linalg.norm(t2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--ham-dir", default="rhf_hamiltonians")
    ap.add_argument("--names", nargs="+", default=["C2H3N_rxn2857_TS", "C2H3N_rxn2858_P"])
    ap.add_argument("--restarts", type=int, default=10)
    ap.add_argument("--kick", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--n-procs", type=int, default=20)
    ap.add_argument("--out", default="gauge_study/exp12_valley_energy.json")
    args = ap.parse_args()
    jobs = [(n, args.data_dir, args.ham_dir, s, args.kick, args.lam, args.maxiter)
            for n in args.names for s in range(-1, args.restarts)]
    cores = sorted(os.sched_getaffinity(0))
    t0 = time.time(); res = []
    with Pool(args.n_procs, initializer=_pin, initargs=(cores,)) as pool:
        for r in pool.imap_unordered(job, jobs):
            res.append(r)
            print(f"  [{len(res)}/{len(jobs)}] {r['name']} seed={r['seed']:3d} resid={r['resid']:.3f} |Z|={r['znorm']:.2f} "
                  f"corr%={100*r['corr_frac']:6.1f} ({r['t']:.0f}s, {time.time()-t0:.0f}s)", flush=True)
    json.dump([{k: v for k, v in r.items() if k != "t2hat"} for r in res], open(args.out, "w"), indent=1)
    print(f"\n=== valley energies (square, lambda={args.lam}, {args.restarts} restarts) ===")
    for n in args.names:
        rs = sorted([r for r in res if r["name"] == n and r["seed"] >= 0], key=lambda r: r["resid"])
        ini = [r for r in res if r["name"] == n and r["seed"] < 0][0]
        cf = 100 * np.array([r["corr_frac"] for r in rs]); rd = np.array([r["resid"] for r in rs])
        T = np.array([r["t2hat"] for r in rs]); D = np.sqrt(((T[:, None] - T[None]) ** 2).sum(-1)) / rs[0]["t2norm"]
        print(f"  {n}: init resid {ini['resid']:.3f} corr% {100*ini['corr_frac']:.1f}")
        print(f"     restarts: resid {rd.min():.3f}-{rd.max():.3f} (median {np.median(rd):.3f}); "
              f"corr% min {cf.min():.1f} median {np.median(cf):.1f} max {cf.max():.1f} (spread {cf.max()-cf.min():.1f} points); "
              f"pairwise dt2hat median {np.median(D[np.triu_indices(len(rs),1)]):.3f}")
        print(f"     corr(resid, corr%) = {np.corrcoef(rd, cf)[0,1]:+.2f}")
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
