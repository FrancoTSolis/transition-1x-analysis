#!/usr/bin/env python3
"""Tables of the TN-reward validation (v2: density-matrix engine), from tn_validate.py jsonl files.

  accuracy  : err = E_TN - E_exact per task and chi; summary per (method, dtype, chi): median / worst error (mHa and
              % of E_HF - E_CCSD), split by norb
  ranking   : per perturbation group and (method, chi): Spearman, Kendall tau, pairwise order agreement, best-member
              agreement, error offset / spread vs the exact spread (noise/signal); aggregated over groups
  n29       : energy vs chi (no exact energy), step per chi doubling, extrapolation in the discarded weight
  timing    : wall / state / <H> per norb and chi; peak GPU memory
Usage: python3 pretrain/rl/tests/tn_v2_report.py --acc a.jsonl [..] --rank r.jsonl [..] --group-json g.json [..]
           [--n29 n.jsonl ..] [--ref-file hi.jsonl --ref-chi 256]  (n29 group self-consistency)
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict

import numpy as np
from scipy.stats import kendalltau, pearsonr, spearmanr


def load(paths):
    out = []
    for p in paths or []:
        for line in open(p):
            if line.strip():
                out.append(json.loads(line))
    return out


def cfg_of(r):
    m = r.get("method", "zipup")
    z = r.get("zip_margin")
    tag = m if m == "dm" else f"zip{z:g}"
    return f"{tag}/{str(r.get('dtype', '')).replace('torch.', '').replace('complex', 'c')}"


def accuracy(recs):
    rows = [r for r in recs if r.get("E_exact") is not None]
    by = defaultdict(dict)
    for r in rows:
        by[(r["norb"], r["name"], r["cand"], round(r["corr_pct_exact"], 2))][(cfg_of(r), r["chi"])] = r
    keys = sorted({k for d in by.values() for k in d}, key=lambda x: (x[0], x[1]))
    print("== accuracy vs exact: err = E_TN - E_exact in mHa (% of E_HF - E_CCSD)")
    print(f"{'norb':4s} {'molecule':18s} {'params':22s} {'exact%':>6s} | " +
          " | ".join(f"{c} chi{x}" for c, x in keys))
    for (norb, name, cand, cex), d in sorted(by.items()):
        cells = []
        for k in keys:
            r = d.get(k)
            cells.append(" " * 14 if r is None else f"{r['err_mHa']:6.2f} ({-r['err_corr_pct']:4.2f})")
        print(f"{norb:4d} {name:18s} {cand[:22]:22s} {cex:6.1f} | " + " | ".join(cells))
    print("\n== accuracy summary (all tasks with exact E)")
    print(f"{'config':10s} {'chi':>4s} {'norb':>6s} {'n':>3s} | {'median err':>10s} {'mean err':>9s} {'worst err':>9s} "
          f"| {'median %':>8s} {'worst %':>7s} | {'>1.6mHa':>7s} | {'max|err| c64-c128?':>5s}")
    for (c, x) in keys:
        for grp, sel in (("all", lambda r: True), ("15-16", lambda r: r["norb"] <= 16), ("17-18", lambda r: r["norb"] >= 17)):
            v = [d[(c, x)] for d in by.values() if (c, x) in d and sel(d[(c, x)])]
            if not v:
                continue
            e = np.array([r["err_mHa"] for r in v])
            p = np.array([-r["err_corr_pct"] for r in v])
            print(f"{c:10s} {x:4d} {grp:>6s} {len(v):3d} | {np.median(e):10.2f} {e.mean():9.2f} {e.max():9.2f} | "
                  f"{np.median(p):8.2f} {p.max():7.2f} | {np.mean(np.abs(e) > 1.6):7.2f} |")


def rank_metrics(E, R):
    E, R = np.asarray(E), np.asarray(R)
    n = len(E)
    err = (E - R) * 1e3
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    agree = np.mean([np.sign(E[i] - E[j]) == np.sign(R[i] - R[j]) for i, j in pairs])
    return {"n": n, "spearman": spearmanr(E, R).correlation, "kendall": kendalltau(E, R).correlation,
            "pearson": pearsonr(E, R)[0], "pair_agree": float(agree), "best_same": bool(np.argmin(E) == np.argmin(R)),
            "ref_std": float(R.std() * 1e3), "err_mean": float(err.mean()), "err_std": float(err.std()),
            "nsr": float(err.std() / max(R.std() * 1e3, 1e-12))}


def ranking(recs, ref_file=None, ref_chi=None, ref_cfg=None):
    ref = {}
    if ref_file:
        for r in load([ref_file]):
            if r["chi"] == ref_chi and (ref_cfg is None or cfg_of(r) == ref_cfg):
                ref[(r["name"], r["cand"])] = r["E"]
    groups = defaultdict(dict)
    for r in recs:
        if "#" not in r["cand"]:
            continue
        key = (r["name"], r["cand"])
        R = ref.get(key) if ref_file else r.get("E_exact")
        if R is None:
            continue
        if ref_file and r["chi"] == ref_chi and (ref_cfg is None or cfg_of(r) == ref_cfg):
            continue
        g = (r["name"], r["cand"].split("#")[0], cfg_of(r), r["chi"])
        groups[g][r["cand"]] = (r["E"], R)
    title = f"reference = chi {ref_chi} ({ref_cfg or 'any'})" if ref_file else "reference = exact (ffsim / GPU exact)"
    print(f"\n== GRPO ranking inside perturbation groups, {title}")
    print(f"{'molecule':18s} {'params':16s} {'config':10s} {'chi':>4s} n | spear kendall pairs best | "
          f"ref std  err mean +- std  noise/signal")
    agg = defaultdict(list)
    for (name, cand, c, x), g in sorted(groups.items(), key=lambda kv: (kv[0][2], kv[0][3], kv[0][0], kv[0][1])):
        keys = sorted(g)
        if len(keys) < 5:
            continue
        m = rank_metrics([g[k][0] for k in keys], [g[k][1] for k in keys])
        agg[(c, x)].append(m)
        print(f"{name:18s} {cand[:16]:16s} {c:10s} {x:4d} {m['n']} | {m['spearman']:5.3f} {m['kendall']:6.3f} "
              f"{m['pair_agree']:5.3f} {str(m['best_same'])[0]:>4s} | {m['ref_std']:6.2f}  {m['err_mean']:+7.2f} +- "
              f"{m['err_std']:5.2f}   {m['nsr']:5.2f}")
    print(f"\n{'config':10s} {'chi':>4s} groups | spearman median (min)  kendall median (min)  pairs mean  best-same  "
          f"noise/signal median (max)")
    for (c, x), ms in sorted(agg.items()):
        sp = np.array([m["spearman"] for m in ms])
        kt = np.array([m["kendall"] for m in ms])
        pa = np.array([m["pair_agree"] for m in ms])
        bs = np.array([m["best_same"] for m in ms])
        ns = np.array([m["nsr"] for m in ms])
        print(f"{c:10s} {x:4d} {len(ms):6d} | {np.median(sp):6.3f} ({sp.min():5.3f})       {np.median(kt):6.3f} "
              f"({kt.min():5.3f})      {pa.mean():6.3f}     {bs.mean():5.2f}      {np.median(ns):5.2f} ({ns.max():4.2f})")


def n29(recs):
    print("\n== norb >= 19: energy vs chi (no exact energy; corr% = (E_HF - E) / (E_HF - E_CCSD))")
    g = defaultdict(dict)
    for r in recs:
        if r["norb"] < 19 or "#" in r["cand"]:
            continue
        g[(r["name"], r["cand"], cfg_of(r))][r["chi"]] = r
    for (name, cand, c), d in sorted(g.items()):
        prev = None
        xs, es, dws = [], [], []
        for x in sorted(d):
            r = d[x]
            step = "" if prev is None else f" step {1e3 * (r['E'] - prev):+7.2f} mHa"
            print(f"{name:18s} {cand:6s} {c:9s} chi {x:4d}: E {r['E']:.6f}  corr {r['corr_pct']:6.2f}%  "
                  f"disc {r['discarded_sum']:.2e}  wall {r['wall']:6.0f}s (state {r['t_state']:5.0f} <H> "
                  f"{r['t_expect']:4.0f})  GPU {r.get('peak_gpu_MB') or 0:6.0f} MB{step}")
            prev = r["E"]
            xs.append(x)
            es.append(r["E"])
            dws.append(r["discarded_sum"])
        if len(xs) >= 3:
            for p in (0.5, 1.0):
                A = np.vstack([np.ones(3), np.array(dws[-3:]) ** p]).T
                coef = np.linalg.lstsq(A, np.array(es[-3:]), rcond=None)[0]
                print(f"    extrapolation E(disc^{p:g} -> 0) from chi {xs[-3:]}: {coef[0]:.6f} "
                      f"(last chi {1e3 * (es[-1] - coef[0]):+.2f} mHa above)")


def timing(recs):
    print("\n== timing (median over records; wall = state + <H>; GPU = peak allocated)")
    g = defaultdict(list)
    for r in recs:
        g[(cfg_of(r), r["norb"], r["chi"])].append(r)
    print(f"{'config':10s} {'norb':>4s} {'chi':>4s} {'n':>4s} | {'wall s':>7s} {'state s':>7s} {'env s':>6s} "
          f"{'sweep s':>7s} {'<H> s':>6s} | {'GPU MB':>7s} | load")
    for (c, nb, x), v in sorted(g.items()):
        med = lambda k: float(np.median([r.get(k) or 0.0 for r in v]))  # noqa: E731
        print(f"{c:10s} {nb:4d} {x:4d} {len(v):4d} | {med('wall'):7.1f} {med('t_state'):7.1f} {med('t_env'):6.1f} "
              f"{med('t_sweep'):7.1f} {med('t_expect'):6.1f} | {med('peak_gpu_MB'):7.0f} | {med('loadavg'):.0f}")


def compare_old(new_recs, old_recs):
    """Old zip-up engine (e.g. the verifier's runs: method missing -> zipup, zip_margin recorded) vs the new default
    on the same (task, chi): error in mHa."""
    print("\n== old zip-up engine vs new default, same task and chi: err mHa (old -> new)")
    new = {}
    for r in new_recs:
        if r.get("E_exact") is not None and "#" not in r["cand"]:
            new[(r["name"], r["cand"], r["chi"])] = r
    rows = defaultdict(dict)
    for r in old_recs:
        if r.get("E_exact") is None or r.get("method", "zipup") != "zipup":
            continue
        key = (r["name"], r["cand"], r["chi"])
        if key in new:
            rows[key][r.get("zip_margin")] = r["err_mHa"]
    for (name, cand, chi), d in sorted(rows.items(), key=lambda kv: (kv[0][2], kv[0][0])):
        olds = "  ".join(f"zip{m:g}: {e:6.2f}" for m, e in sorted(d.items(), key=lambda x: x[0] or 0))
        print(f"{name:18s} {cand[:18]:18s} chi {chi:4d} | {olds:32s} | dm: {new[(name, cand, chi)]['err_mHa']:6.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", nargs="*", default=[], help="old-engine jsonl (zip-up) to compare with --acc")
    ap.add_argument("--acc", nargs="*", default=[])
    ap.add_argument("--rank", nargs="*", default=[])
    ap.add_argument("--n29", nargs="*", default=[])
    ap.add_argument("--ref-file", default=None)
    ap.add_argument("--ref-chi", type=int, default=None)
    ap.add_argument("--ref-cfg", default=None)
    ap.add_argument("--timing", nargs="*", default=[])
    args = ap.parse_args()
    if args.acc:
        accuracy(load(args.acc))
        if args.old:
            compare_old(load(args.acc), load(args.old))
    if args.rank:
        recs = load(args.rank)
        if args.ref_file:
            ranking(recs, args.ref_file, args.ref_chi, args.ref_cfg)
        else:
            ranking(recs)
    if args.n29:
        n29(load(args.n29))
    if args.timing:
        timing(load(args.timing))


if __name__ == "__main__":
    main()
