#!/usr/bin/env python3
"""[Independent verifier] ranking fidelity of TN energies inside GRPO groups.

Reference = exact ffsim energies from verify_tn_pert_refs.py json (--exact-json), or the TN energies of another
(chi, zip_margin, dtype) configuration (--ref-chi / --ref-file) for n29 self-consistency.
Metrics per group and configuration: Spearman, Kendall tau, Pearson (== correlation of GRPO advantages), pairwise
order agreement, sign agreement of the advantages, error offset/std vs exact std, top-1 / bottom-1 agreement.
Usage: python3 pretrain/rl/tests/verify_tn_rank.py --exact-json vpert.json grp_*.jsonl
       python3 pretrain/rl/tests/verify_tn_rank.py --ref-file n29_ref.jsonl --ref-chi 128 n29_*.jsonl
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict

import numpy as np
from scipy.stats import kendalltau, pearsonr, spearmanr


def metrics(E, R):
    E, R = np.asarray(E), np.asarray(R)
    err = (E - R) * 1e3
    n = len(E)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    agree = np.mean([np.sign(E[i] - E[j]) == np.sign(R[i] - R[j]) for i, j in pairs])
    aE, aR = -(E - E.mean()), -(R - R.mean())
    return {"n": n, "spearman": spearmanr(E, R).correlation, "kendall": kendalltau(E, R).correlation,
            "pearson": pearsonr(E, R)[0], "pair_agree": agree, "adv_sign_agree": float(np.mean(np.sign(aE) == np.sign(aR))),
            "best_same": bool(np.argmin(E) == np.argmin(R)), "worst_same": bool(np.argmax(E) == np.argmax(R)),
            "ref_std_mHa": float(R.std() * 1e3), "err_mean_mHa": float(err.mean()), "err_std_mHa": float(err.std()),
            "noise_to_signal": float(err.std() / (R.std() * 1e3))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--exact-json", default=None)
    ap.add_argument("--ref-file", default=None)
    ap.add_argument("--ref-chi", type=int, default=None)
    ap.add_argument("--ref-zip", type=float, default=None)
    args = ap.parse_args()
    recs = [json.loads(l) for f in args.files for l in open(f) if l.strip()]
    ref = {}
    if args.exact_json:
        J = json.load(open(args.exact_json))["per_group"]
        for key, lst in J.items():
            name, cand = key.split("|")
            for k, e in enumerate(lst):
                if e is not None:
                    ref[(name, f"{cand}#pert{k}")] = e["E"]
    if args.ref_file:
        for l in open(args.ref_file):
            r = json.loads(l)
            if r["chi"] == args.ref_chi and (args.ref_zip is None or abs(r["zip_margin"] - args.ref_zip) < 1e-9):
                ref[(r["name"], r["cand"])] = r["E"]
    groups = defaultdict(dict)
    for r in recs:
        key = (r["name"], r["cand"])
        if key not in ref:
            continue
        cfg = (r["name"], r["cand"].split("#")[0], r["chi"], r.get("zip_margin"), r.get("dtype"))
        if args.ref_file and r["chi"] == args.ref_chi and (args.ref_zip is None or abs(r["zip_margin"] - args.ref_zip) < 1e-9):
            continue
        groups[cfg][r["cand"]] = (r["E"], ref[key])
    for cfg, g in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][3] or 0)):
        keys = sorted(g)
        if len(keys) < 3:
            continue
        m = metrics([g[k][0] for k in keys], [g[k][1] for k in keys])
        print(f"{cfg[0]:18s} {cfg[1]:10s} chi {cfg[2]:4d} zip {cfg[3]} {str(cfg[4]).replace('torch.', ''):10s} n={m['n']} "
              f"spearman {m['spearman']:.3f} kendall {m['kendall']:.3f} pearson(=adv corr) {m['pearson']:.3f} "
              f"pairs {m['pair_agree']:.3f} adv-sign {m['adv_sign_agree']:.2f} best {m['best_same']} worst {m['worst_same']} | "
              f"ref std {m['ref_std_mHa']:.2f} mHa, err {m['err_mean_mHa']:+.2f} +- {m['err_std_mHa']:.2f} mHa "
              f"(noise/signal {m['noise_to_signal']:.2f})")


if __name__ == "__main__":
    main()
