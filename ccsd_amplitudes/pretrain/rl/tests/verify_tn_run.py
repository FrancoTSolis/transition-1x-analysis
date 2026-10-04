#!/usr/bin/env python3
"""[Independent verifier] run the PUBLIC API LUCJEnergyTN (default settings unless overridden) on
  * stored tasks with exact references:   --tasks pkl:cand:name,...      (norb 15-17)
  * my perturbation groups:               --group-npz <vpert>.npz [--group-json <vpert>.json] [--members 0,1,..]
  * n29 parameter sets:                   --n29 name:tag,...  (bench_params_<name>.npz in --n29-dir, t1 from rhf_dataset)
  * my own n29 perturbation group:        --n29-group name:tag:seed:G
and append one json line per (task, chi) to --out (E, exact E if known, error, corr%, info, peak GPU memory).
Usage: CUDA_VISIBLE_DEVICES=4 python3 pretrain/rl/tests/verify_tn_run.py --tasks all4:all4_t4:C2H3N_rxn2857_R \
          --chis 32,64 --out <scratch>/acc.jsonl
"""
from __future__ import annotations

import os

_THR = os.environ.get("VTN_THREADS", "2")
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ[_v] = _THR
import argparse  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_energy as T  # noqa: E402
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402

N29_DIR = "/tmp/claude-23917/-xuanwu-tank-east-fts-projects-transition-1x-analysis/f16f5d12-46a8-4d5b-8f5d-48bcb61a4c89/scratchpad"


def exact_refs():
    ref = {}
    for f in ("baselines", "all1", "all4", "all8", "swap"):
        p = ROOT / "pretrain/opt_true/results" / f"energy_small_{f}.json"
        for name, d in json.load(open(p))["per_molecule"].items():
            for cand, v in d.items():
                ref[(name, cand)] = v["E"]
    p = ROOT / "pretrain/opt_true/results/energy_n17_rl1.json"
    for name, d in json.load(open(p))["per_molecule"].items():
        for cand, v in d.items():
            ref[(name, cand)] = v["E"]
    return ref


def load_tasks(args):
    tasks = []
    if args.tasks:
        ref = exact_refs()
        cache = {}
        for spec in args.tasks.split(","):
            pkl, cand, name = spec.split(":")
            if pkl not in cache:
                d = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{pkl}.pkl", "rb"))
                cache[pkl] = {(t[0][0], t[0][1]): t for t in d["tasks"]}
            _, _, U, Z, t1 = cache[pkl][(name, cand)]
            tasks.append((name, cand, U, Z, t1, ref.get((name, cand))))
    if args.group_npz:
        P = np.load(args.group_npz)
        J = json.load(open(args.group_json))["per_group"] if args.group_json else {}
        keys = sorted({k.split("::")[0] for k in P.files})
        for key in keys:
            name, cand = key.split("|")
            if args.group_names and name not in args.group_names.split(","):
                continue
            GU, GZ, t1 = P[f"{key}::U"], P[f"{key}::Z"], P[f"{key}::t1"]
            mem = range(GU.shape[0]) if not args.members else [int(x) for x in args.members.split(",")]
            for k in mem:
                e = J.get(key, [None] * GU.shape[0])[k]
                tasks.append((name, f"{cand}#pert{k}", GU[k], GZ[k], t1, None if e is None else e["E"]))
    if args.n29:
        for spec in args.n29.split(","):
            name, tag = spec.split(":")
            Pp = np.load(Path(args.n29_dir) / f"bench_params_{name}.npz")
            t1 = np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t1"].astype(np.float64)
            tasks.append((name, tag, Pp[f"U_{tag}"], Pp[f"Z_{tag}"], t1, None))
    if args.n29_group:
        from scipy.linalg import expm
        for spec in args.n29_group.split(","):
            name, tag, seed, G = spec.split(":")
            Pp = np.load(Path(args.n29_dir) / f"bench_params_{name}.npz")
            t1 = np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t1"].astype(np.float64)
            U0 = np.stack([T.polar_unitary(u) for u in Pp[f"U_{tag}"].astype(np.complex128)])
            Z0 = Pp[f"Z_{tag}"].astype(np.float64)
            rng = np.random.default_rng([int(seed), 5])
            n = U0.shape[1]
            i = np.arange(n)
            mask = (np.abs(i[:, None] - i[None, :]) <= 1).astype(float)
            mem = [int(x) for x in args.members.split(",")] if args.members else range(int(G) + 1)
            members = [(U0, Z0)]
            for _ in range(int(G)):
                Up, Zp = np.empty_like(U0), Z0.copy()
                for r in range(U0.shape[0]):
                    A = rng.normal(0.0, 0.01, (n, n))
                    Up[r] = U0[r] @ expm(np.triu(A, 1) - np.triu(A, 1).T)
                    B = rng.normal(0.0, 0.005, (n, n))
                    Zp[r] = Z0[r] + (np.triu(B) + np.triu(B, 1).T) * mask
                members.append((Up, Zp))
            for k in mem:
                tasks.append((name, f"{tag}#s{seed}pert{k}", members[k][0], members[k][1], t1, None))
    return tasks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", default=None)
    ap.add_argument("--group-npz", default=None)
    ap.add_argument("--group-json", default=None)
    ap.add_argument("--group-names", default=None)
    ap.add_argument("--members", default=None)
    ap.add_argument("--n29", default=None)
    ap.add_argument("--n29-group", default=None)
    ap.add_argument("--n29-dir", default=N29_DIR)
    ap.add_argument("--chis", default="32,64")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--dtype", default=None, help="complex64/complex128 (default: API default for the device)")
    ap.add_argument("--zip-margin", type=float, default=None)
    ap.add_argument("--cutoff", type=float, default=None)
    ap.add_argument("--basis", default=None, help="None (API default: boys w/ name) | er | pm | boys")
    ap.add_argument("--b2-threads", type=int, default=2)
    ap.add_argument("--basis-cache", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    tasks = load_tasks(args)
    chis = [int(c) for c in args.chis.split(",")]
    scratch_root = Path(args.out).parent / "b2scratch"
    scratch_root.mkdir(parents=True, exist_ok=True)
    evs = {}
    for name, cand, U, Z, t1, E_ex in tasks:
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
        if name not in evs:
            evs.clear()                       # one evaluator alive at a time (see verify_tn_multi.py)
            kw = {}
            if args.dtype:
                kw["dtype"] = getattr(torch, args.dtype)
            if args.zip_margin is not None:
                kw["zip_margin"] = args.zip_margin
            if args.cutoff is not None:
                kw["cutoff"] = args.cutoff
            if args.basis is not None:
                kw["basis"] = args.basis
            t0 = time.time()
            evs[name] = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec,
                                       max_bond=chis[0], device=args.device, name=name,
                                       block2_threads=args.b2_threads, basis_cache=args.basis_cache,
                                       scratch=tempfile.mkdtemp(prefix="vrun_", dir=scratch_root), **kw)
            t_init = time.time() - t0
        ev = evs[name]
        for chi in chis:
            if torch.cuda.is_available():
                torch.cuda.reset_peak_memory_stats()
            t0 = time.time()
            E, info = ev.energy(U, Z, t1, max_bond=chi)
            wall = time.time() - t0
            rec = {"name": name, "cand": cand, "norb": norb, "nelec": list(nelec), "chi": chi, "E": E,
                   "E_exact": E_ex, "err_mHa": None if E_ex is None else (E - E_ex) * 1e3,
                   "corr_pct": (e_hf - E) / (e_hf - e_ccsd) * 100,
                   "corr_pct_exact": None if E_ex is None else (e_hf - E_ex) / (e_hf - e_ccsd) * 100,
                   "err_corr_pct": None if E_ex is None else (E_ex - E) / (e_hf - e_ccsd) * 100,
                   "wall": wall, "t_init": t_init, "device": args.device, "dtype": str(ev.dtype),
                   "zip_margin": ev.zip_margin, "cutoff": ev.cutoff,
                   "peak_gpu_MB": (torch.cuda.max_memory_allocated() / 2 ** 20) if torch.cuda.is_available() else None,
                   **{k: v for k, v in info.items() if not isinstance(v, list)}}
            t_init = 0.0
            print(json.dumps({k: (round(v, 7) if isinstance(v, float) else v) for k, v in rec.items()
                              if k in ("name", "cand", "chi", "E", "E_exact", "err_mHa", "corr_pct", "err_corr_pct",
                                       "wall", "t_state", "t_expect", "t_mpo_build", "discarded_sum", "max_bond",
                                       "peak_gpu_MB", "imag")}), flush=True)
            with open(args.out, "a") as f:
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
