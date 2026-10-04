#!/usr/bin/env python3
"""TN energy vs exact on stored reference tasks (norb 15-17): error (Hartree, % corr) vs max_bond, timings.

Usage: CUDA_VISIBLE_DEVICES=4 python3 pretrain/rl/tests/tn_bench_small.py --pkl all1 --cand all1_t1 \
          --names C2H3N_rxn2857_P --chis 64,128,256 --out pretrain/rl/tests/results/bench_small.jsonl
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ.setdefault(_v, "8")
import argparse  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import re  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_energy as T  # noqa: E402
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402


def exact_refs():
    """(name, cand) -> exact E from the stored energy jobs (small: json; norb 17: logs + json)."""
    ref = {}
    for f in ("baselines", "all1", "all4", "all8", "swap"):
        p = ROOT / "pretrain/opt_true/results" / f"energy_small_{f}.json"
        if p.exists():
            for name, d in json.load(open(p))["per_molecule"].items():
                for cand, v in d.items():
                    ref[(name, cand)] = v["E"]
    p = ROOT / "pretrain/opt_true/results/energy_n17_rl1.json"
    if p.exists():
        for name, d in json.load(open(p))["per_molecule"].items():
            for cand, v in d.items():
                ref[(name, cand)] = v["E"]
    return ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", default="all1")
    ap.add_argument("--cand", default="all1_t1")
    ap.add_argument("--names", required=True)
    ap.add_argument("--chis", default="64,128,256")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--dtype", default="complex64")
    ap.add_argument("--b2-threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--perturbed", default=None, help="npz from tn_perturbed_refs.py (uses its params instead)")
    args = ap.parse_args()
    dtype = getattr(torch, args.dtype)
    ref = exact_refs()
    tasks = []
    if args.perturbed:
        P = np.load(args.perturbed, allow_pickle=True)
        pj = json.load(open(args.perturbed.replace(".npz", ".json")))["per_molecule"]
        for i, name in enumerate([str(x) for x in P["names"]]):
            if name not in args.names.split(","):
                continue
            for k in range(P["U"].shape[1]):
                E = pj[name][k]["E"] if pj[name][k] is not None else None
                tasks.append((name, f"pert{k}", P["U"][i, k], P["Z"][i, k], P["t1"][i], E))
    else:
        d = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{args.pkl}.pkl", "rb"))
        by = {(t[0][0], t[0][1]): t for t in d["tasks"]}
        for name in args.names.split(","):
            _, _, U, Z, t1 = by[(name, args.cand)]
            tasks.append((name, args.cand, U, Z, t1, ref.get((name, args.cand))))
    chis = [int(c) for c in args.chis.split(",")]
    ev_cache = {}
    for name, cand, U, Z, t1, E_ex in tasks:
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
        if name not in ev_cache:
            ev_cache.clear()
            ev_cache[name] = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec,
                                            device=args.device, dtype=dtype, block2_threads=args.b2_threads)
        ev = ev_cache[name]
        for chi in chis:
            E, info = ev.energy(U, Z, t1, max_bond=chi)
            rec = {"name": name, "cand": cand, "norb": norb, "chi": chi, "E": E, "E_exact": E_ex,
                   "err_mHa": None if E_ex is None else (E - E_ex) * 1e3,
                   "corr_frac": (e_hf - E) / (e_hf - e_ccsd),
                   "corr_frac_exact": None if E_ex is None else (e_hf - E_ex) / (e_hf - e_ccsd),
                   "err_corr_pct": None if E_ex is None else (E - E_ex) / (e_hf - e_ccsd) * -100,
                   "dtype": args.dtype, **{k: v for k, v in info.items() if k != "bond_dims"}}
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)
            with open(args.out, "a") as f:
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
