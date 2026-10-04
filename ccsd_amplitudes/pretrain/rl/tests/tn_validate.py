#!/usr/bin/env python3
"""TN reward validation runner: LUCJEnergyTN on tasks with (or without) exact references; one json line per
(task, chi) with EVERY engine setting recorded (method, zip_margin, cutoff, dtype, device, basis, mode_tol, ...).

Task sources (combinable):
  --tasks pkl:cand:name,...      stored energy tasks (runs_ot/energy_tasks/<pkl>.pkl); exact E from the stored result
                                 JSONs (pretrain/opt_true/results: energy_small_*, energy_n17_rl1, energy_largeval_*,
                                 energy_smallval_rl4L)
  --groups f.npz [--group-json f.json] [--group-keys k1,k2] [--members 0,1,..]
                                 perturbation groups: tn_group_refs.py / verify_tn_pert_refs.py format
                                 ("<name>|<cand>::U") or tn_perturbed_refs.py format (names, U (G, M, ...))
  --n29 name:tag,...             n29 parameter sets (<n29-dir>/bench_params_<name>.npz, t1 from rhf_dataset)
  --n29-group name:tag:seed:G    n29 perturbation group (generator of verify_tn_run.py)
Usage: CUDA_VISIBLE_DEVICES=4 python3 pretrain/rl/tests/tn_validate.py --tasks all1:all1_t1:C2H3N_rxn2857_P \
          --chis 64,128 --out pretrain/rl/tests/results/v2_acc.jsonl
"""
from __future__ import annotations

import os

_THR = os.environ.get("TNV_THREADS", "2")
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ[_v] = _THR
import argparse  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import socket  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import importlib  # noqa: E402

T = importlib.import_module(os.environ.get("TN_MODULE", "pretrain.rl.tn_energy"))
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402
from pretrain.rl.tests.tn_group_refs import stored_exact  # noqa: E402

N29_DIR = "/tmp/claude-23917/-xuanwu-tank-east-fts-projects-transition-1x-analysis/f16f5d12-46a8-4d5b-8f5d-48bcb61a4c89/scratchpad"


def load_groups(path, json_path, keys, members):
    P = np.load(path, allow_pickle=True)
    out = []
    if "names" in P.files:                                     # tn_perturbed_refs.py format
        J = json.load(open(json_path or path.replace(".npz", ".json")))["per_molecule"]
        cand = str(P["cand"]) if "cand" in P.files else "c"
        for i, name in enumerate([str(x) for x in P["names"]]):
            key = f"{name}|{cand}"
            if keys and key not in keys and name not in keys:
                continue
            for k in range(P["U"].shape[1]):
                if members and k not in members:
                    continue
                e = J[name][k]
                out.append((name, f"{cand}#pert{k}", P["U"][i, k], P["Z"][i, k], P["t1"][i],
                            None if e is None else e["E"]))
        return out
    J = json.load(open(json_path))["per_group"] if json_path else {}
    for key in sorted({k.split("::")[0] for k in P.files}):
        name, cand = key.split("|")
        if keys and key not in keys and name not in keys:
            continue
        GU, GZ, t1 = P[f"{key}::U"], P[f"{key}::Z"], P[f"{key}::t1"]
        t1 = None if t1.size == 0 else t1
        for k in range(GU.shape[0]):
            if members and k not in members:
                continue
            e = J.get(key, [None] * GU.shape[0])[k]
            out.append((name, f"{cand}#pert{k}", GU[k], GZ[k], t1, None if e is None else e["E"]))
    return out


def load_tasks(args):
    tasks = []
    members = [int(x) for x in args.members.split(",")] if args.members else None
    if args.tasks:
        ref = stored_exact()
        cache = {}
        for spec in args.tasks.split(","):
            pkl, cand, name = spec.split(":")
            if pkl not in cache:
                d = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{pkl}.pkl", "rb"))
                cache[pkl] = {(t[0][0], t[0][1]): t for t in d["tasks"]}
            _, _, U, Z, t1 = cache[pkl][(name, cand)]
            tasks.append((name, cand, U, Z, t1, ref.get((name, cand))))
    for g in (args.groups or []):
        keys = args.group_keys.split(",") if args.group_keys else None
        tasks += load_groups(g, args.group_json, keys, members)
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
            mem = [(U0, Z0)]
            for _ in range(int(G)):
                Up, Zp = np.empty_like(U0), Z0.copy()
                for r in range(U0.shape[0]):
                    A = rng.normal(0.0, 0.01, (n, n))
                    Up[r] = U0[r] @ expm(np.triu(A, 1) - np.triu(A, 1).T)
                    B = rng.normal(0.0, 0.005, (n, n))
                    Zp[r] = Z0[r] + (np.triu(B) + np.triu(B, 1).T) * mask
                mem.append((Up, Zp))
            for k in (members if members else range(int(G) + 1)):
                tasks.append((name, f"{tag}#s{seed}pert{k}", mem[k][0], mem[k][1], t1, None))
    return tasks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", default=None)
    ap.add_argument("--groups", action="append", default=None)
    ap.add_argument("--group-json", default=None)
    ap.add_argument("--group-keys", default=None)
    ap.add_argument("--members", default=None)
    ap.add_argument("--n29", default=None)
    ap.add_argument("--n29-group", default=None)
    ap.add_argument("--n29-dir", default=N29_DIR)
    ap.add_argument("--chis", default="64")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--dtype", default=None, help="complex64 / complex128 (default: API default)")
    ap.add_argument("--method", default=None, help="dm / zipup (default: API default)")
    ap.add_argument("--zip-margin", type=float, default=None)
    ap.add_argument("--cutoff", type=float, default=None)
    ap.add_argument("--basis", default=None)
    ap.add_argument("--b2-threads", type=int, default=1)
    ap.add_argument("--stack-mem-gb", type=float, default=None)
    ap.add_argument("--basis-cache", default=str(ROOT / "pretrain/rl/tests/results/split_bases"))
    ap.add_argument("--repeat", type=int, default=1, help="evaluate every (task, chi) this many times (timing)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    tasks = load_tasks(args)
    chis = [int(c) for c in args.chis.split(",")]
    scratch_root = Path(tempfile.gettempdir()) / f"tnv_{os.getpid()}"
    scratch_root.mkdir(parents=True, exist_ok=True)
    gpu = torch.cuda.get_device_name(0) if (args.device != "cpu" and torch.cuda.is_available()) else None
    cur, ev = None, None
    for name, cand, U, Z, t1, E_ex in tasks:
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
        if cur != name:
            ev = None
            kw = {}
            for k_, v_ in (("method", args.method), ("zip_margin", args.zip_margin), ("cutoff", args.cutoff),
                           ("basis", args.basis)):
                if v_ is not None:
                    kw[k_] = v_
            if args.dtype:
                kw["dtype"] = getattr(torch, args.dtype)
            if args.stack_mem_gb is not None:
                kw["stack_mem"] = int(args.stack_mem_gb * (1 << 30))
            t0 = time.time()
            ev = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec,
                                max_bond=chis[0], device=args.device, name=name, block2_threads=args.b2_threads,
                                basis_cache=args.basis_cache,
                                scratch=tempfile.mkdtemp(prefix="b2_", dir=scratch_root), **kw)
            t_init = time.time() - t0
            cur = name
        for chi in chis:
            for rep in range(args.repeat):
                if torch.cuda.is_available():
                    torch.cuda.reset_peak_memory_stats()
                load0 = os.getloadavg()[0]
                t0 = time.time()
                E, info = ev.energy(U, Z, t1, max_bond=chi)
                wall = time.time() - t0
                rec = {"name": name, "cand": cand, "norb": norb, "nelec": list(nelec), "chi": chi, "rep": rep,
                       "E": E, "E_exact": E_ex, "err_mHa": None if E_ex is None else (E - E_ex) * 1e3,
                       "corr_pct": (e_hf - E) / (e_hf - e_ccsd) * 100,
                       "corr_pct_exact": None if E_ex is None else (e_hf - E_ex) / (e_hf - e_ccsd) * 100,
                       "err_corr_pct": None if E_ex is None else (E_ex - E) / (e_hf - e_ccsd) * 100,
                       "e_hf": e_hf, "e_ccsd": e_ccsd, "wall": wall, "t_init": t_init,
                       "peak_gpu_MB": (torch.cuda.max_memory_allocated() / 2 ** 20) if torch.cuda.is_available() else None,
                       "host": socket.gethostname(), "gpu": gpu, "loadavg": load0, "threads": _THR,
                       **ev.settings(), **{k: v for k, v in info.items() if not isinstance(v, (list, dict))}}
                t_init = 0.0
                print(json.dumps({k: (round(v, 7) if isinstance(v, float) else v) for k, v in rec.items()
                                  if k in ("name", "cand", "chi", "method", "zip_margin", "dtype", "E", "err_mHa",
                                           "corr_pct", "wall", "t_state", "t_expect", "discarded_sum",
                                           "discarded_zip_sum", "max_bond", "peak_gpu_MB")}), flush=True)
                with open(args.out, "a") as f:
                    f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
