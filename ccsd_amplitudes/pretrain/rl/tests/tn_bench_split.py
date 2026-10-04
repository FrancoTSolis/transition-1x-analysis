#!/usr/bin/env python3
"""Split-localized-basis TN energy (LUCJEnergySplitTN) vs exact references: error vs max_bond, timings.

Tasks: stored energy-task pickles (--pkl/--cand/--names) or perturbation groups (--perturbed npz from
tn_perturbed_refs.py), or n29 parameter sets (--n29-params dir with bench_params_<name>.npz, --tags net,frame).
Usage: CUDA_VISIBLE_DEVICES=4 python3 pretrain/rl/tests/tn_bench_split.py --names C2H3N_rxn2857_P --chis 64,128 \
          --out pretrain/rl/tests/results/split_bench.jsonl
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ.setdefault(_v, "6")
import argparse  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_energy as T  # noqa: E402
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402
from pretrain.rl.tests.tn_bench_small import exact_refs  # noqa: E402


def load_tasks(args):
    tasks = []
    names = args.names.split(",")
    if args.perturbed:
        P = np.load(args.perturbed, allow_pickle=True)
        pj = json.load(open(args.perturbed.replace(".npz", ".json")))["per_molecule"]
        for i, name in enumerate([str(x) for x in P["names"]]):
            if names != ["all"] and name not in names:
                continue
            for k in range(P["U"].shape[1]):
                if args.members and k not in [int(x) for x in args.members.split(",")]:
                    continue
                E = pj[name][k]["E"] if pj[name][k] is not None else None
                tasks.append((name, f"pert{k}", P["U"][i, k], P["Z"][i, k], P["t1"][i], E))
    elif args.n29_params:
        for name in names:
            Pp = np.load(Path(args.n29_params) / f"bench_params_{name}.npz")
            t1 = np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t1"].astype(np.float64)
            for tag in args.tags.split(","):
                tasks.append((name, tag, Pp[f"U_{tag}"], Pp[f"Z_{tag}"], t1, None))
    else:
        ref = exact_refs()
        d = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{args.pkl}.pkl", "rb"))
        by = {(t[0][0], t[0][1]): t for t in d["tasks"]}
        for name in names:
            for cand in args.cand.split(","):
                _, _, U, Z, t1 = by[(name, cand)]
                tasks.append((name, cand, U, Z, t1, ref.get((name, cand))))
    return tasks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", default="all1")
    ap.add_argument("--cand", default="all1_t1")
    ap.add_argument("--names", required=True)
    ap.add_argument("--perturbed", default=None)
    ap.add_argument("--members", default=None, help="comma list of group member indices (default: all)")
    ap.add_argument("--n29-params", default=None)
    ap.add_argument("--tags", default="net")
    ap.add_argument("--chis", default="64,128")
    ap.add_argument("--method", default="boys")
    ap.add_argument("--order", default="fiedler")
    ap.add_argument("--cutoff", type=float, default=1e-12)
    ap.add_argument("--mode-tol", type=float, default=1e-8)
    ap.add_argument("--zip-margin", type=float, default=2.0)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--dtype", default="complex64")
    ap.add_argument("--b2-threads", type=int, default=6)
    ap.add_argument("--basis-cache", default=str(ROOT / "pretrain/rl/tests/results/split_bases"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    dtype = getattr(torch, args.dtype)
    tasks = load_tasks(args)
    chis = [int(c) for c in args.chis.split(",")]
    cur, ev = None, None
    for name, cand, U, Z, t1, E_ex in tasks:
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
        if cur != name:
            t0 = time.time()
            S, occ = T.split_localized_basis(name, ROOT, args.method, args.order, cache_dir=args.basis_cache)
            t_basis = time.time() - t0
            ev = T.LUCJEnergySplitTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec, S, occ,
                                     cutoff=args.cutoff, mode_tol=args.mode_tol, device=args.device, dtype=dtype,
                                     block2_threads=args.b2_threads)
            ev.zip_margin = args.zip_margin
            cur = name
        for chi in chis:
            E, info = ev.energy(U, Z, t1, max_bond=chi)
            rec = {"name": name, "cand": cand, "norb": norb, "chi": chi, "E": E, "E_exact": E_ex,
                   "err_mHa": None if E_ex is None else (E - E_ex) * 1e3,
                   "corr_frac": (e_hf - E) / (e_hf - e_ccsd),
                   "corr_frac_exact": None if E_ex is None else (e_hf - E_ex) / (e_hf - e_ccsd),
                   "err_corr_pct": None if E_ex is None else (E_ex - E) / (e_hf - e_ccsd) * 100,
                   "method": args.method, "order": args.order, "dtype": args.dtype, "t_basis": t_basis,
                   "cutoff": args.cutoff, "mode_tol": args.mode_tol, **info}
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)
            with open(args.out, "a") as f:
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
