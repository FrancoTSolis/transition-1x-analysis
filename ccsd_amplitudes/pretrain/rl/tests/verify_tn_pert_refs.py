#!/usr/bin/env python3
"""[Independent verifier] exact ffsim energies of NEW GRPO-like perturbation groups (own generator, own seeds).

Group around a stored task (pkl, cand, name): member 0 = polar(U), Z (centre; must reproduce the stored exact E);
members k >= 1:  U'_r = polar(U_r) @ expm(K_r), K_r real antisymmetric, upper entries iid N(0, sig_k^2);
                 Z'_r = Z_r + E_r, E_r symmetric, entries iid N(0, sig_z^2) on |i - j| <= 1, 0 elsewhere.
One thread per worker (CPU rules).  Writes <out>.npz and <out>.json (incrementally).
Usage: python3 pretrain/rl/tests/verify_tn_pert_refs.py --spec all1:label:C2H3N_rxn2858_P,all4:all4_t4:C2H3N_rxn2857_R \
          --group 8 --seed 777 --n-procs 7 --out <scratch>/vpert
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
from scipy.linalg import expm  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))


def polar(U):
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def make_group(U, Z, rng, G, sig_k, sig_z):
    n_reps, n, _ = U.shape
    i = np.arange(n)
    mask = (np.abs(i[:, None] - i[None, :]) <= 1).astype(float)
    Us, Zs = [U.copy()], [Z.copy()]
    for _ in range(G):
        Up = np.empty_like(U)
        Zp = Z.copy()
        for r in range(n_reps):
            A = rng.normal(0.0, sig_k, (n, n))
            K = np.triu(A, 1) - np.triu(A, 1).T
            Up[r] = U[r] @ expm(K)
            B = rng.normal(0.0, sig_z, (n, n))
            Esym = np.triu(B) + np.triu(B, 1).T
            Zp[r] = Z[r] + Esym * mask
        Us.append(Up)
        Zs.append(Zp)
    return np.stack(Us), np.stack(Zs)


def _work(a):
    key, k, name, U, Z, t1 = a
    from pretrain.rl.energy import exact_energy, make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
    t0 = time.time()
    E = exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1=t1))
    return key, k, E, (e_hf - E) / (e_hf - e_ccsd), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True, help="comma list pkl:cand:name")
    ap.add_argument("--group", type=int, default=8)
    ap.add_argument("--seed", type=int, default=777)
    ap.add_argument("--sig-k", type=float, default=0.01)
    ap.add_argument("--sig-z", type=float, default=0.005)
    ap.add_argument("--n-procs", type=int, default=7)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    specs = [s.split(":") for s in args.spec.split(",")]
    jobs, groups = [], {}
    for gi, (pkl, cand, name) in enumerate(specs):
        d = pickle.load(open(ROOT / "runs_ot" / "energy_tasks" / f"{pkl}.pkl", "rb"))
        by = {(t[0][0], t[0][1]): t for t in d["tasks"]}
        _, _, U, Z, t1 = by[(name, cand)]
        U = np.stack([polar(u) for u in np.asarray(U, dtype=np.complex128)])
        Z = np.asarray(Z, dtype=np.float64)
        t1 = np.asarray(t1, dtype=np.float64)
        rng = np.random.default_rng([args.seed, gi, 99])
        GU, GZ = make_group(U, Z, rng, args.group, args.sig_k, args.sig_z)
        key = f"{name}|{cand}"
        groups[key] = {"name": name, "cand": cand, "pkl": pkl, "U": GU, "Z": GZ, "t1": t1}
        jobs += [(key, k, name, GU[k], GZ[k], t1) for k in range(args.group + 1)]
    np.savez(args.out + ".npz", **{f"{key}::{f}": g[f] for key, g in groups.items() for f in ("U", "Z", "t1")})
    res = {key: [None] * (args.group + 1) for key in groups}
    t0 = time.time()
    # interleave so that the centre members finish first (convention check early)
    jobs.sort(key=lambda j: (j[1], j[0]))
    with Pool(args.n_procs) as pool:
        for key, k, E, cf, dt in pool.imap_unordered(_work, jobs):
            res[key][k] = {"E": E, "corr_frac": cf, "t": dt}
            print(f"  {key} k={k} E={E:.10f} corr%={cf * 100:.4f} ({dt:.0f}s, {time.time() - t0:.0f}s total)",
                  flush=True)
            json.dump({"args": vars(args), "per_group": res}, open(args.out + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
