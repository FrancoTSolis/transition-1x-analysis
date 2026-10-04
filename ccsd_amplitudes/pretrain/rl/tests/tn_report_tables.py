#!/usr/bin/env python3
"""Tables for the TN-reward feasibility report, from pretrain/rl/tests/results/*.jsonl."""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

R = Path(__file__).resolve().parent / "results"
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


def load(*names):
    out = []
    for n in names:
        p = R / n
        if p.exists():
            out += [json.loads(l) for l in open(p) if l.strip()]
    return out


def main():
    print("== accuracy vs exact (norb 15-17); err = E_TN - E_exact; CPU runs complex128, GPU runs complex64")
    rows = load("split_margin_1.5.jsonl", "split_n15_gpu.jsonl", "split_acc_n15_16.jsonl", "split_acc_n17.jsonl")
    seen = set()
    by = defaultdict(dict)
    for r in rows:
        key = (r["name"], r["cand"], r["chi"])
        if key in seen:
            continue
        seen.add(key)
        by[(r["name"], r["cand"], r["norb"], round(100 * r["corr_frac_exact"], 2))][r["chi"]] = r
    chis = sorted({c for d in by.values() for c in d})
    print(f"{'molecule':18s} {'params':8s} norb corr_ex% | " + " | ".join(f"chi {c:3d}: err mHa (%corr) t_s" for c in chis))
    for (name, cand, norb, cex), d in sorted(by.items(), key=lambda x: (x[0][2], x[0][0])):
        cells = []
        for c in chis:
            if c in d:
                r = d[c]
                cells.append(f"{r['err_mHa']:7.2f} ({-r['err_corr_pct']:5.2f}) {r['t_state'] + r['t_expect']:5.0f}")
            else:
                cells.append(" " * 27)
        print(f"{name:18s} {cand:8s} {norb:4d} {cex:8.2f} | " + " | ".join(cells))
    print("\n== GRPO ranking (exact references, norb 15, chi 32)")
    from pretrain.rl.tests.tn_rank_eval import main as rank
    rank([str(R / "split_pert_np_chi32.jsonl")])
    print("\n== norb 29 energy vs chi (no exact; corr% = (E_HF - E)/(E_HF - E_CCSD))")
    n29 = load("split_n29_gpu.jsonl")
    g = defaultdict(dict)
    for r in n29:
        g[(r["name"], r["cand"])][r["chi"]] = r
    for (name, cand), d in g.items():
        prev = None
        for c in sorted(d):
            r = d[c]
            step = "" if prev is None else f"  step {1e3 * (r['E'] - prev):+.2f} mHa"
            print(f"{name:18s} {cand:5s} chi {c:4d}: E {r['E']:.6f}  corr {100 * r['corr_frac']:6.2f}%  disc {r['discarded_sum']:.1e}"
                  f"  t_state {r['t_state']:5.0f}s  t_expect {r['t_expect']:5.1f}s  maxbond {r['max_bond']}{step}")
            prev = r["E"]
    if (R / "split_pert_n29_chi64.jsonl").exists():
        print("\n== norb 29 GRPO ranking self-consistency (chi 32 vs chi 64 as reference)")
        rank([str(R / "split_pert_n29_chi32.jsonl")], ref=str(R / "split_pert_n29_chi64.jsonl"))


if __name__ == "__main__":
    main()
