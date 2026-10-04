#!/usr/bin/env python3
"""Exact (ffsim statevector) energies of GRPO-like perturbation groups, for the TN-reward ranking test.

For each molecule, the centre is a stored task (runs_ot/energy_tasks/<pkl>, candidate <cand>); U is re-unitarized
(SVD polar factor) as in the energy jobs.  Group members k = 1..G:
    U'_r = U_r @ expm(K_r),  K_r real antisymmetric, independent entries ~ N(0, sig_k^2)   (per rep r)
    Z'_r = Z_r + E_r,        E_r symmetric ~ N(0, sig_z^2) on the square-mask entries (diagonal + first off-diagonal)
Member 0 is the unperturbed centre (its energy is compared with the stored exact value as a convention check).

Writes <out>.npz (names, U (M,G+1,2,n,n), Z, t1) and <out>.json (per molecule list of E, corr_frac, time).
CPU rules: one thread per worker (OMP/MKL/OPENBLAS/RAYON pinned before numpy import), --n-procs workers.

Usage: python3 -m pretrain.rl.tests.tn_perturbed_refs --names C2H3N_rxn2857_P,C2H3N_rxn2857_TS --pkl all1 \
           --cand all1_t1 --group 8 --n-procs 8 --out pretrain/rl/tests/results/perturbed_n15
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))


def polar(U):
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def square_mask(n):
    i = np.arange(n)
    return np.abs(i[:, None] - i[None, :]) <= 1


def perturb(U, Z, rng, sig_k=0.01, sig_z=0.005):
    from scipy.linalg import expm
    n_reps, n, _ = U.shape
    Up, Zp = np.empty_like(U), Z.copy()
    m = square_mask(n)
    for r in range(n_reps):
        K = np.triu(rng.normal(0.0, sig_k, (n, n)), 1)
        K = K - K.T
        Up[r] = U[r] @ expm(K)
        E = np.triu(rng.normal(0.0, sig_z, (n, n)))
        E = E + np.triu(E, 1).T
        Zp[r] = Z[r] + E * m
    return Up, Zp


def _work(args):
    name, k, U, Z, t1 = args
    from pretrain.rl.energy import exact_energy, make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
    t0 = time.time()
    E = exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1=t1))
    return name, k, E, (e_hf - E) / (e_hf - e_ccsd), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True)
    ap.add_argument("--pkl", default="all1")
    ap.add_argument("--cand", default="all1_t1")
    ap.add_argument("--group", type=int, default=8)
    ap.add_argument("--sig-k", type=float, default=0.01)
    ap.add_argument("--sig-z", type=float, default=0.005)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--n-procs", type=int, default=8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    names = args.names.split(",")
    d = pickle.load(open(ROOT / "runs_ot" / "energy_tasks" / f"{args.pkl}.pkl", "rb"))
    by = {(t[0][0], t[0][1]): t for t in d["tasks"]}
    Us, Zs, T1s = [], [], []
    jobs = []
    for i, name in enumerate(names):
        _, _, U, Z, t1 = by[(name, args.cand)]
        U = polar(np.asarray(U, dtype=np.complex128))
        Z = np.asarray(Z, dtype=np.float64)
        rng = np.random.default_rng([args.seed, i])
        grpU, grpZ = [U], [Z]
        for _ in range(args.group):
            Up, Zp = perturb(U, Z, rng, args.sig_k, args.sig_z)
            grpU.append(Up)
            grpZ.append(Zp)
        Us.append(np.stack(grpU))
        Zs.append(np.stack(grpZ))
        T1s.append(np.asarray(t1, dtype=np.float64))
        jobs += [(name, k, grpU[k], grpZ[k], T1s[-1]) for k in range(args.group + 1)]
    np.savez(args.out + ".npz", names=np.array(names), U=np.stack(Us), Z=np.stack(Zs),
             t1=np.stack(T1s) if len({t.shape for t in T1s}) == 1 else np.array(T1s, dtype=object),
             cand=args.cand, pkl=args.pkl, sig_k=args.sig_k, sig_z=args.sig_z, seed=args.seed)
    res = {n: [None] * (args.group + 1) for n in names}
    t0 = time.time()
    with Pool(args.n_procs) as pool:
        for name, k, E, cf, dt in pool.imap_unordered(_work, jobs):
            res[name][k] = {"E": E, "corr_frac": cf, "t": dt}
            print(f"  {name} k={k} E={E:.8f} corr%={cf * 100:.3f} ({dt:.0f}s, {time.time() - t0:.0f}s total)",
                  flush=True)
            json.dump({"args": vars(args), "per_molecule": res}, open(args.out + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
