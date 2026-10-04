#!/usr/bin/env python3
"""Tables of docs/rl_larger_molecules.md from the result files (prints markdown).

  python3 -m pretrain.rl.report_tables
Missing result files are skipped.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain" / "opt_true" / "results"


def load(name):
    p = RES / name
    return json.load(open(p))["per_molecule"] if p.exists() else None


def merge(*pms):
    out = {}
    for pm in pms:
        for n, d in (pm or {}).items():
            out.setdefault(n, {}).update({c: v["corr_frac"] for c, v in d.items() if not v.get("error")})
    return out


def log_evals(run, kind="eval_val"):
    p = ROOT / "rl_runs" / run / "log.jsonl"
    if not p.exists():
        return {}
    rows = [json.loads(ln) for ln in open(p)]
    return {r["step"]: r for r in rows if r.get("type") == kind}


def table(per, cands, labels, ref="label", names=None):
    names = names or sorted(per)
    names = [n for n in names if ref in per[n]]
    r = np.array([per[n][ref] for n in names])
    print(f"| candidate | n | mean % corr | median | min | vs {ref} (mean) | better than {ref} |")
    print("|:--|--:|--:|--:|--:|--:|--:|")
    for c in cands:
        ok = [i for i, n in enumerate(names) if c in per[n]]
        if not ok:
            continue
        v = np.array([per[names[i]][c] for i in ok])
        d = v - r[ok]
        cmp_ = "—" if c == ref else f"{100 * d.mean():+.1f}"
        bt = "—" if c == ref else f"{int((d > 0).sum())}/{len(ok)}"
        print(f"| {labels.get(c, c)} | {len(ok)} | {100 * v.mean():.1f} | {100 * np.median(v):.1f} | "
              f"{100 * v.min():.1f} | {cmp_} | {bt} |")
    print()


def curve(run, steps=None):
    ev, tr = log_evals(run), log_evals(run, "eval_train")
    if not ev:
        return
    steps = steps or sorted(ev)
    print("| step | " + " | ".join(str(s) for s in steps) + " |")
    print("|:--|" + "--:|" * len(steps))
    print("| val mean % corr | " + " | ".join(f"{100 * ev[s]['corr_mean']:.1f}" for s in steps) + " |")
    print("| val median | " + " | ".join(f"{100 * ev[s]['corr_median']:.1f}" for s in steps) + " |")
    print("| train mean | " + " | ".join(f"{100 * tr[s]['corr_mean']:.1f}" if s in tr else "" for s in steps) + " |")
    print()


LAB = {"label": "optimize=True label", "pre1": "pretrained, one shot", "pre4": "pretrained, 4 recycles",
       "rl1": "small-set RL, one shot", "rl1s": "small-set RL, one shot", "rl4": "small-set RL, 4 recycles",
       "rl4s": "small-set RL, 4 recycles", "rl4L": "norb 15-18 RL (step 25)", "rl4L_s50": "norb 15-18 RL (step 50)",
       "rl4Lf": "norb 15-18 RL (step 50)", "rl4n29": "+ n29 TN RL (best chi-64 val step)",
       "rl4n29f": "+ n29 TN RL (step 30)"}


def convergence(n29):
    """Same molecules at three accuracies: does the candidate ranking depend on the truncation?"""
    c128 = merge(load("energy_n29val_labrl4L_chi128.json"), load("energy_n29val_rl4n29_chi128.json"))
    c3 = merge(load("energy_n29val_sub6_chi256m3.json"), load("energy_n29val_rl4n29_sub6_chi256m3.json"))
    if not (n29 and c128 and c3):
        return
    names = sorted(n for n in c3 if "label" in c3[n] and "rl4L" in c3[n])
    cands = [c for c in ["label", "rl4L", "rl4n29"] if all(c in c3[n] and c in c128.get(n, {}) for n in names)]
    print(f"## n29 truncation check ({len(names)} val molecules, mean % corr)\n")
    print("| MPS | " + " | ".join(LAB.get(c, c) for c in cands) + " | norb 15-18 RL − label |")
    print("|:--|" + "--:|" * (len(cands) + 1))
    for title, pm in [("chi 128, margin 1.5", c128), ("chi 256, margin 1.5", n29), ("chi 256, margin 3", c3)]:
        v = {c: np.array([pm[n][c] for n in names]) for c in cands}
        print(f"| {title} | " + " | ".join(f"{100 * v[c].mean():.2f}" for c in cands)
              + f" | {100 * (v['rl4L'] - v['label']).mean():+.2f} |")
    print()


def main():
    print("## norb 15-18 RL curve (rl_runs/grpo_large_rl4)\n")
    curve("grpo_large_rl4")
    per = merge(load("energy_largeval_base.json"), load("energy_largeval_pre4_rl4L.json"),
                load("energy_largeval_rl4n29.json"), load("energy_largeval_rl4n29_C2H6O.json"))
    ev = log_evals("grpo_large_rl4")
    if per and 50 in ev:
        for n, v in ev[50]["per_mol"].items():
            per.setdefault(n, {})["rl4L_s50"] = v
    if per:
        print("## norb 17-18 held-out (22 molecules, exact)\n")
        table(per, ["label", "pre1", "pre4", "rl1", "rl4", "rl4L", "rl4L_s50", "rl4n29", "rl4n29f"], LAB)
    sm = merge(load("energy_smallval_rl4L.json"), load("energy_smallval_rl4n29.json"))
    if sm:
        print("## original small val (9 molecules, norb 15-16, exact)\n")
        table(sm, ["rl4s", "rl4L", "rl4n29", "rl4n29f"], LAB, ref="rl4s")
    n19 = merge(load("energy_n19_all.json"), load("energy_n19_rl4n29.json"))
    if n19:
        print("## norb 19 test (12 C2H3NO, exact; never used in RL or selection)\n")
        table(n19, ["label", "pre1", "pre4", "rl1s", "rl4s", "rl4L", "rl4Lf", "rl4n29", "rl4n29f"], LAB)
    n29 = merge(load("energy_n29val_base_chi256.json"), load("energy_n29val_rl4n29_chi256.json"))
    if n29:
        print("## n29 val (20 molecules, MPS chi 256, zip margin 1.5)\n")
        table(n29, ["label", "pre1", "pre4", "rl4s", "rl4L", "rl4n29", "rl4n29f"], LAB)
    for f, title in [("energy_n29val_labrl4L_chi128.json", "chi 128"),
                     ("energy_n29val_sub6_chi256m3.json", "chi 256, zip margin 3 (6-molecule subset)")]:
        pm = merge(load(f))
        if pm:
            print(f"## n29 val at {title}\n")
            table(pm, ["label", "rl4L", "rl4n29"], LAB)
    convergence(n29)
    te = merge(load("energy_n29test_base_chi256.json"), load("energy_n29test_rl4n29_chi256.json"))
    if te:
        print("## n29 test (24 molecules from 24 unused reactions, MPS chi 256)\n")
        table(te, ["label", "pre4", "rl4s", "rl4L", "rl4n29", "rl4n29f"], LAB)
        print("By formula (mean % corr; C4H5NO and C3H5N3 are the n29 RL formulas):\n")
        forms = sorted({n.split("_")[0] for n in te})
        cands = [c for c in ["label", "pre4", "rl4s", "rl4L", "rl4n29", "rl4n29f"] if any(c in te[n] for n in te)]
        print("| formula (nocc, nvirt) | n | " + " | ".join(LAB.get(c, c) for c in cands) + " |")
        print("|:--|--:|" + "--:|" * len(cands))
        idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
        for fm in forms:
            ns = [n for n in te if n.split("_")[0] == fm]
            sh = tuple(idx[ns[0]][1:])
            print(f"| {fm} {sh} | {len(ns)} | " + " | ".join(
                f"{100 * np.mean([te[n][c] for n in ns if c in te[n]]):.1f}" if any(c in te[n] for n in ns) else ""
                for c in cands) + " |")
        print()
    print("## n29 TN RL curve (rl_runs/grpo_n29_tn, chi 64 rewards and evals)\n")
    curve("grpo_n29_tn")


if __name__ == "__main__":
    main()
