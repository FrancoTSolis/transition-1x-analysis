#!/usr/bin/env python3
"""Experiment 10: how many distinct compressed-DF minima does a molecule have?

For the multiple-hypothesis (winner-takes-all) amortized optimizer we need to
know how many basins matter.  Per molecule: M restarts of ffsim's compressed
optimization (square, lambda=0.005, 500 it) from the canonical init with
random column phases and a small random anti-Hermitian kick of scale --kick on
kappa.  Solutions are clustered by the relative distance of their invariant
content t2hat (threshold --thr); we report #clusters, their sizes, the
residual of each cluster's best member, and the residual spread.

Run from ccsd_amplitudes/ (ffsim venv):
    python3 -m gauge_study.exp10_basins --n-mols 4 --restarts 10 --n-procs 10
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
    name, data_dir, seed, kick, lam, maxiter = args
    d = np.load(Path(data_dir) / f"{name}.npz")
    t2 = d["t2"].astype(np.float64)
    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    init = canonical_exact_init(t2)
    m = z_mask("square", norb)[None]
    rng = np.random.default_rng(seed)
    U0 = init.U.copy()
    if seed > 0:
        ph = np.exp(2j * np.pi * rng.random((2, norb)))
        U0 = U0 * ph[:, None, :]
        for k in range(2):
            A = rng.standard_normal((norb, norb)) + 1j * rng.standard_normal((norb, norb))
            K = kick * (A - A.conj().T) / np.sqrt(2 * norb)
            U0[k] = U0[k] @ scipy.linalg.expm(K)
    r = compress_from_init(t2, init.Z * m, U0, connectivity="square", regularization=lam,
                           znorm_ref=init.znorm_full, maxiter=maxiter)
    th = reconstruct_t2(r["Z"], r["U"], nocc)
    return dict(name=name, seed=seed, resid=rel_residual(t2, r["Z"], r["U"]),
                znorm=float(np.linalg.norm(r["Z"])), t2hat=th, t2norm=float(np.linalg.norm(t2)),
                resid_init=rel_residual(t2, init.Z * m, init.U))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--names-file", default="gauge_study/names_n29_16_13.txt")
    ap.add_argument("--n-mols", type=int, default=4)
    ap.add_argument("--restarts", type=int, default=10)
    ap.add_argument("--kick", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--thr", type=float, default=0.05, help="cluster threshold on ||dt2hat||/||t2||")
    ap.add_argument("--n-procs", type=int, default=10)
    ap.add_argument("--out", default="gauge_study/exp10_basins.json")
    args = ap.parse_args()
    names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    rng = np.random.default_rng(1)
    names = [names[i] for i in rng.choice(len(names), args.n_mols, replace=False)]
    jobs = [(n, args.data_dir, s, args.kick, args.lam, args.maxiter) for n in names for s in range(args.restarts)]
    cores = sorted(os.sched_getaffinity(0))
    t0 = time.time()
    res = []
    with Pool(args.n_procs, initializer=_pin, initargs=(cores,)) as pool:
        for r in pool.imap_unordered(job, jobs):
            res.append(r)
            print(f"  [{len(res)}/{len(jobs)}] {r['name']} seed={r['seed']} resid={r['resid']:.3f} |Z|={r['znorm']:.2f} ({time.time()-t0:.0f}s)", flush=True)
    summary = {}
    print(f"\n=== basins per molecule (square, lambda={args.lam}, {args.restarts} restarts, kick={args.kick}, thr={args.thr}) ===")
    for n in names:
        rs = sorted([r for r in res if r["name"] == n], key=lambda r: r["resid"])
        T = np.array([r["t2hat"].ravel() for r in rs]); nt = rs[0]["t2norm"]
        D = np.sqrt(((T[:, None] - T[None]) ** 2).sum(-1)) / nt
        # greedy clustering by t2hat distance
        labels = -np.ones(len(rs), int); c = 0
        for i in range(len(rs)):
            if labels[i] < 0:
                labels[(D[i] < args.thr) & (labels < 0)] = c; c += 1
        sizes = [int((labels == k).sum()) for k in range(c)]
        best = [min(r["resid"] for r, l in zip(rs, labels) if l == k) for k in range(c)]
        resids = [r["resid"] for r in rs]
        print(f"  {n:24s} init {rs[0]['resid_init']:.3f} | resid min {min(resids):.3f} median {np.median(resids):.3f} max {max(resids):.3f}"
              f" | clusters {c} sizes {sizes} best-resid {['%.3f' % b for b in best]} | pairwise dt2hat median {np.median(D[np.triu_indices(len(rs),1)]):.3f}")
        summary[n] = dict(resids=resids, clusters=c, sizes=sizes, best=best,
                          dt2hat_median=float(np.median(D[np.triu_indices(len(rs), 1)])), resid_init=rs[0]["resid_init"])
    json.dump(summary, open(args.out, "w"), indent=1)
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
