#!/usr/bin/env python3
"""Ranking fidelity of an approximate energy inside GRPO groups.

Input: jsonl records with name, cand (pert0..pertG), chi, E, E_exact (tn_bench_split.py --perturbed ...).
Per molecule and chi: Spearman and Pearson correlation between approximate and exact energies over the group,
the spread of the exact energies (std, mHa), the error offset (mean) and the error spread (std, mHa) -- the part of
the error that can reorder group members -- and the correlation of GRPO advantages (r - mean)/std with r = -E.
Usage: python3 pretrain/rl/tests/tn_rank_eval.py results/split_pert_chi32.jsonl [more.jsonl ...]
       python3 pretrain/rl/tests/tn_rank_eval.py --ref results/split_pert_n29_chi64.jsonl results/split_pert_n29_chi32.jsonl
       (--ref: no exact energies (norb 29); the reference file's E plays the role of E_exact -> self-consistency)
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import pearsonr, spearmanr


def main(paths, ref=None):
    recs = [json.loads(l) for p in paths for l in open(p) if l.strip()]
    if ref is not None:
        rmap = {}
        for l in open(ref):
            if l.strip():
                r = json.loads(l)
                rmap[(r["name"], r["cand"])] = (r["E"], r["chi"])
        for r in recs:
            if (r["name"], r["cand"]) in rmap:
                r["E_exact"] = rmap[(r["name"], r["cand"])][0]
                r["variant"] = f"vs chi{rmap[(r['name'], r['cand'])][1]}"
    groups = defaultdict(dict)
    for r in recs:
        if r.get("E_exact") is None:
            continue
        groups[(r["name"], r["chi"], r.get("variant", ""))][r["cand"]] = r
    out = []
    for (name, chi, var), g in sorted(groups.items()):
        keys = sorted(g)
        E = np.array([g[k]["E"] for k in keys])
        Ex = np.array([g[k]["E_exact"] for k in keys])
        err = (E - Ex) * 1e3
        rho = spearmanr(E, Ex).correlation if len(keys) > 2 else float("nan")
        pr = pearsonr(E, Ex)[0] if len(keys) > 2 else float("nan")
        adv = lambda x: (x - x.mean()) / (x.std() + 1e-12)  # noqa: E731
        ac = float(np.corrcoef(adv(-E), adv(-Ex))[0, 1]) if len(keys) > 2 else float("nan")
        row = {"name": name, "chi": chi, "variant": var, "n": len(keys), "spearman": rho, "pearson": pr,
               "adv_corr": ac, "exact_std_mHa": float(Ex.std() * 1e3), "err_mean_mHa": float(err.mean()),
               "err_std_mHa": float(err.std()), "err_max_abs_mHa": float(np.abs(err).max())}
        out.append(row)
        print(f"{name:20s} chi {chi:5d} {var:10s} n={len(keys)}  spearman {rho:6.3f}  pearson {pr:6.3f}  adv-corr {ac:6.3f} "
              f"| exact std {row['exact_std_mHa']:.2f} mHa | err mean {row['err_mean_mHa']:+.3f} std "
              f"{row['err_std_mHa']:.3f} max|.| {row['err_max_abs_mHa']:.3f} mHa")
    return out


if __name__ == "__main__":
    args = sys.argv[1:]
    ref = None
    if args and args[0] == "--ref":
        ref, args = args[1], args[2:]
    main(args, ref)
