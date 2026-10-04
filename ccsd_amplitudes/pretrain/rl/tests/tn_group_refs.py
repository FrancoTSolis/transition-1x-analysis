#!/usr/bin/env python3
"""Exact energies of GRPO-like perturbation groups for the TN ranking test (GPU exact engine, pretrain.rl.gpu_energy).

Group around a stored task (pkl:cand:name from runs_ot/energy_tasks/<pkl>.pkl): member 0 = polar(U), Z (the centre;
its energy must reproduce the stored ffsim energy, checked and printed); members k >= 1:
    U'_r = polar(U_r) @ expm(K_r),  K_r real antisymmetric, upper entries iid N(0, sig_k^2)
    Z'_r = Z_r + E_r,  E_r symmetric, entries iid N(0, sig_z^2) on |i - j| <= 1 (the square mask), 0 elsewhere
(the generator of pretrain/rl/tests/verify_tn_pert_refs.py; seeds differ per group: rng([seed, group index, 99])).
Exact energies: LUCJEnergyGPU (complex64 default: |dE| <~ 1e-7 Ha against ffsim, see
pretrain/rl/tests/results/gpu_energy_refs_small_c16.json; --dtype complex128 for an fp64 run).  norb <= 17 on a
12 GB GPU.  Writes <out>.npz (keys "<name>|<cand>::U/Z/t1") and <out>.json (per_group[key] = [{E, corr_frac, t}]).

Usage: CUDA_VISIBLE_DEVICES=4 python3 pretrain/rl/tests/tn_group_refs.py \
           --spec all1:all1_t1:C2H4O_rxn0721_P,swap:all1_t1_swap:C2N2_rxn3923_P --group 8 --seed 2026 \
           --out pretrain/rl/tests/results/groups_n16
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import argparse  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
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
        Up, Zp = np.empty_like(U), Z.copy()
        for r in range(n_reps):
            A = rng.normal(0.0, sig_k, (n, n))
            Up[r] = U[r] @ expm(np.triu(A, 1) - np.triu(A, 1).T)
            B = rng.normal(0.0, sig_z, (n, n))
            Zp[r] = Z[r] + (np.triu(B) + np.triu(B, 1).T) * mask
        Us.append(Up)
        Zs.append(Zp)
    return np.stack(Us), np.stack(Zs)


def stored_exact():
    """(name, cand) -> exact ffsim E from the stored result JSONs (norb 15-18)."""
    ref = {}
    R = ROOT / "pretrain/opt_true/results"
    for f in ["energy_small_baselines.json", "energy_small_all1.json", "energy_small_all4.json",
              "energy_small_all8.json", "energy_small_swap.json", "energy_n17_rl1.json",
              "energy_largeval_base.json", "energy_smallval_rl4L.json", "energy_largeval_pre4_rl4L.json"]:
        p = R / f
        if not p.exists():
            continue
        for name, d in json.load(open(p))["per_molecule"].items():
            for cand, v in d.items():
                if isinstance(v, dict) and v.get("E") is not None and np.isfinite(v["E"]):
                    ref[(name, cand)] = float(v["E"])
    return ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True, help="comma list pkl:cand:name")
    ap.add_argument("--group", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--sig-k", type=float, default=0.01)
    ap.add_argument("--sig-z", type=float, default=0.005)
    ap.add_argument("--dtype", default="complex64")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import torch
    from pretrain.rl.gpu_energy import LUCJEnergyGPU
    from pretrain.rl.hamiltonian import load_hamiltonian
    ref = stored_exact()
    groups = {}
    for gi, spec in enumerate(args.spec.split(",")):
        pkl, cand, name = spec.split(":")
        d = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{pkl}.pkl", "rb"))
        by = {(t[0][0], t[0][1]): t for t in d["tasks"]}
        _, _, U, Z, t1 = by[(name, cand)]
        U = np.stack([polar(u) for u in np.asarray(U, dtype=np.complex128)])
        Z = np.asarray(Z, dtype=np.float64)
        t1 = None if t1 is None else np.asarray(t1, dtype=np.float64)
        rng = np.random.default_rng([args.seed, gi, 99])
        GU, GZ = make_group(U, Z, rng, args.group, args.sig_k, args.sig_z)
        groups[f"{name}|{cand}"] = {"name": name, "cand": cand, "pkl": pkl, "U": GU, "Z": GZ, "t1": t1}
    np.savez(args.out + ".npz", **{f"{k}::{f}": (g[f] if g[f] is not None else np.zeros(0))
                                   for k, g in groups.items() for f in ("U", "Z", "t1")})
    res = {k: [None] * (args.group + 1) for k in groups}
    meta = {}
    dtype = getattr(torch, args.dtype)
    t00 = time.time()
    for key, g in groups.items():
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", g["name"])
        eng = LUCJEnergyGPU(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec, dtype=dtype)
        for k in range(args.group + 1):
            t0 = time.time()
            E = float(eng.energy(g["U"][k], g["Z"][k], t1=g["t1"]))
            res[key][k] = {"E": E, "corr_frac": (e_hf - E) / (e_hf - e_ccsd), "t": time.time() - t0}
        E0 = ref.get((g["name"], g["cand"]))
        Es = np.array([x["E"] for x in res[key]])
        meta[key] = {"norb": norb, "nelec": list(nelec), "e_hf": e_hf, "e_ccsd": e_ccsd, "stored_exact_centre": E0,
                     "centre_minus_stored": None if E0 is None else res[key][0]["E"] - E0,
                     "std_mHa": float(Es.std() * 1e3), "range_mHa": float(np.ptp(Es) * 1e3)}
        print(f"{key:40s} norb {norb} centre {res[key][0]['E']:.9f} stored {E0} "
              f"(diff {meta[key]['centre_minus_stored']}) group std {meta[key]['std_mHa']:.3f} mHa "
              f"({time.time() - t00:.0f}s)", flush=True)
        del eng
        torch.cuda.empty_cache()
        json.dump({"args": vars(args), "engine": f"LUCJEnergyGPU {args.dtype}", "meta": meta, "per_group": res},
                  open(args.out + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
