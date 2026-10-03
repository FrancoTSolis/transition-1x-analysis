#!/usr/bin/env python3
"""Merge exact-energy result files (eval_slot_energy JSON) into one table over the 30 held-out small molecules.

Usage: python3 -m pretrain.opt_true.energy_table pretrain/opt_true/results/energy_small_*.json [--md out.md]
"""
from __future__ import annotations

import argparse
import json

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--md", default=None)
    args = ap.parse_args()
    per, resid = {}, {}
    for f in args.files:
        d = json.load(open(f))
        for n, rows in d["per_molecule"].items():
            for k, v in rows.items():
                per.setdefault(k, {})[n] = 100 * v["corr_frac"]
                resid.setdefault(k, {})[n] = d["meta"][n]["resid"].get(k, float("nan"))
    names = sorted(set().union(*[set(v) for v in per.values()]))
    lab = per.get("label", {})
    order = ["init", "frame", "frame+adam300", "label", "label_swap"] + sorted(k for k in per if k not in
                                                                                ("init", "frame", "frame+adam300", "label", "label_swap"))
    lines = ["| candidate | n | median % corr | mean | min | median residual | vs label: median diff | better than label |",
             "|:--|--:|--:|--:|--:|--:|--:|--:|"]
    for k in order:
        if k not in per:
            continue
        ns = [n for n in names if n in per[k]]
        cf = np.array([per[k][n] for n in ns])
        rs = np.array([resid[k][n] for n in ns])
        both = [n for n in ns if n in lab]
        diff = np.array([per[k][n] - lab[n] for n in both]) if k != "label" and both else np.array([])
        dtxt = f"{np.median(diff):+.1f}" if len(diff) else "—"
        btxt = f"{int((diff > 0).sum())}/{len(diff)}" if len(diff) else "—"
        lines.append(f"| {k} | {len(ns)} | {np.median(cf):.1f} | {cf.mean():.1f} | {cf.min():.1f} | {np.nanmedian(rs):.3f} | {dtxt} | {btxt} |")
    # best of each pair (layer order chosen by energy)
    for base in ("label", "frame+adam300", "all1_t1", "all1_t1+adam100"):
        if base in per and base + "_swap" in per:
            ns = [n for n in names if n in per[base] and n in per[base + "_swap"]]
            best = np.array([max(per[base][n], per[base + "_swap"][n]) for n in ns])
            diff = np.array([max(per[base][n], per[base + "_swap"][n]) - lab[n] for n in ns if n in lab])
            lines.append(f"| best layer order of {base} | {len(ns)} | {np.median(best):.1f} | {best.mean():.1f} | {best.min():.1f} | — | "
                         f"{np.median(diff):+.1f} | {int((diff > 0).sum())}/{len(diff)} |")
    txt = "\n".join(lines)
    print(txt)
    if args.md:
        open(args.md, "w").write(txt + "\n")


if __name__ == "__main__":
    main()
