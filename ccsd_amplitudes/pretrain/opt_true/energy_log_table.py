#!/usr/bin/env python3
"""Per-candidate energy summary from eval_slot_energy *logs* (for jobs stopped before writing their JSON).

Parses lines "  <name>  <candidate>  corr%  <x>  resid <r>  (<t>s)" from one or more logs and compares every
candidate with 'label' on the molecules where both are available.

Usage: python3 -m pretrain.opt_true.energy_log_table rl_runs/energy17_*.out [--names gauge_study/names_norb17_c2hno.txt]
"""
from __future__ import annotations

import argparse
import re

import numpy as np

PAT = re.compile(r"^\s+(\S+)\s+(\S+)\s+corr%\s+(-?[\d.]+)\s+resid\s+([\d.]+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--names", default=None)
    args = ap.parse_args()
    keep = None
    if args.names:
        keep = {ln.strip() for ln in open(args.names) if ln.strip()}
    per, res = {}, {}
    for f in args.logs:
        for ln in open(f):
            m = PAT.match(ln)
            if m:
                n, k, c, r = m.group(1), m.group(2), float(m.group(3)), float(m.group(4))
                if keep is None or n in keep:
                    per.setdefault(k, {})[n] = c
                    res.setdefault(k, {})[n] = r
    lab = per.get("label", {})
    print("| candidate | n | mean % corr | median | min | median residual | paired vs label (mean diff, better) |")
    print("|:--|--:|--:|--:|--:|--:|:--|")
    for k in sorted(per, key=lambda x: (x != "label", x)):
        ns = sorted(per[k])
        c = np.array([per[k][n] for n in ns])
        r = np.array([res[k][n] for n in ns])
        both = [n for n in ns if n in lab]
        if k != "label" and both:
            d = np.array([per[k][n] - lab[n] for n in both])
            pv = f"{d.mean():+.1f} on {len(both)}, better on {int((d > 0).sum())}"
        else:
            pv = "—"
        print(f"| {k} | {len(ns)} | {c.mean():.1f} | {np.median(c):.1f} | {c.min():.1f} | {np.median(r):.3f} | {pv} |")


if __name__ == "__main__":
    main()
