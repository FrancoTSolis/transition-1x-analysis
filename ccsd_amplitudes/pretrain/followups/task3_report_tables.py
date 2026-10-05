#!/usr/bin/env python3
"""Compact headline tables of task 3 (stdout) from the summary.json written by task3_analyze.py.

    python3 pretrain/followups/task3_report_tables.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
S = json.load(open(ROOT / "pretrain/opt_true/results/followups/task3_converged_labels/summary.json"))
SETS = ["norb 15-16 val (9)", "norb 17-18 val (22)", "norb 19 test (12)", "n29 val (20)", "n29 test (24)"]
NET = ["pre4", "rl4L", "rl4n29", "rl4n29f"]


def wl(p):
    return f"{p['wins']}/{p['n']}"


def main():
    print("| set (energies) | label, 500 it. | label, converged | pretrained ×4 | norb 15-18 RL | + n29 RL step 15 "
          "| + n29 RL step 30 |")
    print("|:--|--:|--:|--:|--:|--:|--:|")
    for s in SETS + ["pooled"]:
        if s not in S:
            continue
        t = S[s]["table"] if s != "pooled" else S[s]
        kind = "exact" if s != "pooled" and S[s]["kind"] == "exact" else ("MPS χ256" if s != "pooled" else "mixed")
        n = S[s]["n"] if s != "pooled" else t["label"]["n"]
        cells = [f"{t['label']['mean']:.1f}", f"{t['labconv']['mean']:.1f}"]
        cells += [f"{t[c]['mean']:.1f}" if c in t else "" for c in NET]
        print(f"| {s.split(' (')[0] if s != 'pooled' else 'all sets pooled'} ({n}, {kind}) | " + " | ".join(cells) + " |")
    print()
    print("Paired differences in points of % corr (wins = molecules where the candidate is higher):\n")
    print("| set | converged − 500-it. label | pretrained ×4 − converged | norb 15-18 RL − converged "
          "| step-15 − converged | step-30 − converged |")
    print("|:--|--:|--:|--:|--:|--:|")
    for s in SETS + ["pooled"]:
        if s not in S:
            continue
        t = S[s]["table"] if s != "pooled" else S[s]
        c0 = t["labconv"]["vs_label"]
        cells = [f"{c0['mean']:+.2f} ({c0['wins']} up, {c0['losses']} down)"]
        cells += [f"{t[c]['vs_labconv']['mean']:+.2f} ({wl(t[c]['vs_labconv'])})" if c in t else "" for c in NET]
        print(f"| {s.split(' (')[0] if s != 'pooled' else 'all sets pooled'} | " + " | ".join(cells) + " |")
    print()


if __name__ == "__main__":
    main()
