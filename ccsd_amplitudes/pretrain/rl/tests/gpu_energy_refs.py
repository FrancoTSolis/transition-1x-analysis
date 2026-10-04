#!/usr/bin/env python3
"""Validate pretrain/rl/gpu_energy.py against stored exact (ffsim, CPU) LUCJ energies.

norb 15-16:  runs_ot/energy_tasks/<tag>.pkl  +  pretrain/opt_true/results/energy_small_<tag>.json
norb 17:     runs_ot/energy_tasks/n17_rl1.pkl + pretrain/opt_true/results/energy_n17_rl1.json (full precision)
             runs_ot/energy_tasks/n17_pre.pkl + log lines "<name> <cand> corr% <x>" (1 decimal) in
             rl_runs/energy17_*.out

  CUDA_VISIBLE_DEVICES=2 python -m pretrain.rl.tests.gpu_energy_refs --set small --dtypes c8 c16 \
      --per-mol 2 --out pretrain/rl/tests/results/gpu_energy_refs_small.json
"""
from __future__ import annotations

import argparse
import json
import pickle
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402

from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac  # noqa: E402

SMALL_TAGS = ["baselines", "all1", "all4", "all8", "swap"]


def load_small(tags):
    """[(name, cand, U, Z, t1, E_ref, cf_ref, src)]"""
    out = []
    for tag in tags:
        tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{tag}.pkl", "rb"))["tasks"]
        ref = json.load(open(ROOT / "pretrain/opt_true/results" / f"energy_small_{tag}.json"))["per_molecule"]
        for (name, cand), _, U, Z, t1 in tasks:
            r = ref.get(name, {}).get(cand)
            if r is None:
                continue
            out.append((name, cand, U, Z, t1, float(r["E"]), float(r["corr_frac"]), f"energy_small_{tag}.json"))
    return out


def load_n17():
    out = []
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks/n17_rl1.pkl", "rb"))["tasks"]
    ref = json.load(open(ROOT / "pretrain/opt_true/results/energy_n17_rl1.json"))["per_molecule"]
    for (name, cand), _, U, Z, t1 in tasks:
        r = ref[name][cand]
        out.append((name, cand, U, Z, t1, float(r["E"]), float(r["corr_frac"]), "energy_n17_rl1.json"))
    # 1-decimal log references
    logref = {}
    for f in ["rl_runs/energy17_54608122.out", "rl_runs/energy17_54613781.out"]:
        for ln in open(ROOT / f):
            m = re.match(r"\s+(\S+)\s+(\S+)\s+corr%\s+([-\d.]+)", ln)
            if m:
                logref[(m.group(1), m.group(2))] = (float(m.group(3)) / 100, f)
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks/n17_pre.pkl", "rb"))["tasks"]
    for (name, cand), _, U, Z, t1 in tasks:
        if (name, cand) in logref:
            cf, f = logref[(name, cand)]
            out.append((name, cand, U, Z, t1, float("nan"), cf, Path(f).name + " (1 decimal)"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["small", "n17"], default="small")
    ap.add_argument("--tags", nargs="+", default=SMALL_TAGS)
    ap.add_argument("--dtypes", nargs="+", default=["c8", "c16"])
    ap.add_argument("--per-mol", type=int, default=0, help="max tasks per molecule (0 = all)")
    ap.add_argument("--names", nargs="*", default=None)
    ap.add_argument("--max-mem-gb", type=float, default=None)
    ap.add_argument("--backend", default="auto")
    ap.add_argument("--reuse-engine", action="store_true",
                    help="switch molecules of the same (norb, nelec) with set_npz() instead of a new engine")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = load_small(args.tags) if args.set == "small" else load_n17()
    if args.names:
        rows = [r for r in rows if r[0] in args.names]
    by_mol: dict = {}
    for r in rows:
        by_mol.setdefault(r[0], []).append(r)
    if args.per_mol:
        # rotate through candidate types so the selection mixes them
        sel = {}
        for i, (name, rs) in enumerate(sorted(by_mol.items())):
            rs = sorted(rs, key=lambda r: (r[1], r[7]))
            off = (i * args.per_mol) % len(rs)
            sel[name] = [rs[(off + j) % len(rs)] for j in range(min(args.per_mol, len(rs)))]
        by_mol = sel
    print(f"{sum(len(v) for v in by_mol.values())} tasks on {len(by_mol)} molecules", flush=True)

    results = []
    cache = {}
    for name, rs in sorted(by_mol.items()):
        for dt_s in args.dtypes:
            dt = torch.complex64 if dt_s == "c8" else torch.complex128
            torch.cuda.reset_peak_memory_stats()
            t0 = time.time()
            d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
            key = (int(d["norb"]), int(d["nelec_a"]), dt_s)
            if args.reuse_engine and key in cache:
                eng = cache[key].set_npz(name)
            else:
                for e in cache.values():
                    e.release()
                cache.clear()
                eng = LUCJEnergyGPU.from_npz(name, dtype=dt, max_mem_gb=args.max_mem_gb, backend=args.backend)
                if args.reuse_engine:
                    cache[key] = eng
            t_setup = time.time() - t0
            for (nm, cand, U, Z, t1, E_ref, cf_ref, src) in rs:
                E = eng.energy(U, Z, t1=t1, connectivity="square")
                cf = corr_frac(E, eng.e_hf, eng.e_ccsd)
                res = dict(name=nm, cand=cand, dtype=dt_s, norb=eng.norb, k=eng.k, E=E, E_ref=E_ref,
                           dE=E - E_ref, cf=cf, cf_ref=cf_ref, dcf=cf - cf_ref, src=src,
                           t=eng.timing["total"], timing=eng.timing, norm=eng.last["norm"],
                           t_setup=t_setup, peak_gb=torch.cuda.max_memory_allocated() / 1e9,
                           tile=eng.T, nq=eng.nq)
                results.append(res)
                print(f"  {nm:20s} {cand:24s} {dt_s:3s} n={eng.norb} k={eng.k}  E {E:.10f}  dE {E - E_ref:+.2e}"
                      f"  corr% {100 * cf:8.4f} (ref {100 * cf_ref:8.4f}, d {cf - cf_ref:+.1e})"
                      f"  |psi|^2-1 {eng.last['norm'] - 1:+.1e}  {eng.timing['total']:.2f}s", flush=True)
            if not args.reuse_engine:
                eng.release()
            del eng
            torch.cuda.empty_cache()

    print("\n=== summary ===")
    summ = {}
    for dt_s in args.dtypes:
        rr = [r for r in results if r["dtype"] == dt_s]
        full = [r for r in rr if np.isfinite(r["E_ref"])]
        s = dict(n=len(rr), n_full_precision=len(full))
        if full:
            s.update(max_abs_dE=max(abs(r["dE"]) for r in full), max_abs_dcf=max(abs(r["dcf"]) for r in full),
                     mean_abs_dE=float(np.mean([abs(r["dE"]) for r in full])))
        low = [r for r in rr if not np.isfinite(r["E_ref"])]
        if low:
            s["max_abs_dcf_vs_1decimal"] = max(abs(r["dcf"]) for r in low)
        s["mean_t"] = float(np.mean([r["t"] for r in rr]))
        summ[dt_s] = s
        print(f"  {dt_s}: {json.dumps(s)}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(summary=summ, results=results, args=vars(args)), open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
