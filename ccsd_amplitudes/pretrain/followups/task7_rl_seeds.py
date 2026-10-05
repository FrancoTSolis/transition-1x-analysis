#!/usr/bin/env python3
"""Task 7 (Oct 2026 follow-ups): seed-to-seed variation of the norb 15-18 exact-reward GRPO run.

Seed 0 is the run of docs/rl_larger_molecules.md section 3 (rl_runs/grpo_large_rl4).  Seeds 1, 2, ... repeat its
recipe with only --seed changed (pretrain/followups/task7_launch_seed.sh) and live in
rl_runs/followups/task7_rl_seeds/seed<k>/.

Held-out energies of each seed's best (val-selected) and final (step 50) policy:
  norb 19 (gauge_study/names_norb19.txt, 12, exact)   energy_n19_seed<k>.json      tags s<k>best / s<k>last
  norb 15-16 small val (pretrain/rl/small_val.txt, 9)  energy_smallval_seed<k>.json same tags
Seed 0 uses the existing files (energy_n19_all.json rl4L / rl4Lf, energy_smallval_rl4L.json rl4L) plus
energy_smallval_seed0.json (its step-50 policy on the small val set, tag s0last).

  python3 -m pretrain.followups.task7_rl_seeds eval --seed 1 --kind n19 [--which best last]   # dump + queue energies
  python3 -m pretrain.followups.task7_rl_seeds report [--seeds 0 1 2 3] [--json OUT]   # markdown tables + JSON
"eval" runs the unchanged core tools (pretrain.rl.policy_dump on CPU, pretrain.rl.queue_eval through the shared
exact queue).  When a seed's best checkpoint is its final one, the policy is evaluated once and the result is
stored under both tags.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain" / "opt_true" / "results"
FRES = RES / "followups" / "task7_rl_seeds"
RUNS = ROOT / "rl_runs" / "followups" / "task7_rl_seeds"


def run_dir(seed):
    return ROOT / "rl_runs" / "grpo_large_rl4" if seed == 0 else RUNS / f"seed{seed}"


def read_log(seed):
    p = run_dir(seed) / "log.jsonl"
    if not p.exists():
        return None
    rows = [json.loads(ln) for ln in open(p) if ln.strip()]
    return {"val": {r["step"]: r for r in rows if r.get("type") == "eval_val"},
            "train_eval": {r["step"]: r for r in rows if r.get("type") == "eval_train"},
            "train": {r["step"]: r for r in rows if r.get("type") == "train"},
            "done": any(r.get("type") == "done" for r in rows)}


def best_step(val):
    """The driver's choice: step 0 is the initial best, a later eval replaces it only if strictly better."""
    steps = sorted(val)
    b = steps[0]
    for s in steps[1:]:
        if val[s]["corr_mean"] > val[b]["corr_mean"]:
            b = s
    return b


def per_mol(path, key):
    p = RES / path
    if not p.exists():
        return None
    pm = json.load(open(p))["per_molecule"]
    out = {n: d[key]["corr_frac"] for n, d in pm.items() if key in d and not d[key].get("error")}
    return out or None


def e_hf(name):
    return float(np.load(ROOT / "rhf_dataset" / f"{name}.npz")["e_hf"])


def distinct(names):
    """One name per distinct molecule (the norb-17 and norb-19 sets repeat a reactant; compare by E_HF)."""
    seen, out = set(), []
    for n in names:
        k = round(e_hf(n), 7)
        if k not in seen:
            seen.add(k)
            out.append(n)
    return out


def held_out(seed, which, kind):
    """kind 'n19' | 'small'; which 'best' | 'last' -> {name: corr frac} or None."""
    if seed == 0:
        if kind == "n19":
            return per_mol("energy_n19_all.json", "rl4L" if which == "best" else "rl4Lf")
        if which == "best":
            return per_mol("energy_smallval_rl4L.json", "rl4L")
        return per_mol("followups/task7_rl_seeds/energy_smallval_seed0.json", "s0last")
    f = f"followups/task7_rl_seeds/energy_{'n19' if kind == 'n19' else 'smallval'}_seed{seed}.json"
    return per_mol(f, f"s{seed}{which}")


def stats(v, ref, names):
    """mean / median / min of v and its gain over ref (percentage points) on names."""
    a = np.array([v[n] for n in names])
    d = a - np.array([ref[n] for n in names])
    return {"n": len(names), "mean": 100 * a.mean(), "median": 100 * np.median(a), "min": 100 * a.min(),
            "gain": 100 * d.mean(), "gain_se_mols": 100 * d.std(ddof=1) / np.sqrt(len(d)), "better": int((d > 0).sum())}


def mstd(x):
    x = np.asarray(x, float)
    return (x.mean(), x.std(ddof=1) if len(x) > 1 else float("nan"))


def tci(x, conf=0.95):
    """t confidence interval of the mean across seeds."""
    from scipy import stats as st
    x = np.asarray(x, float)
    if len(x) < 2:
        return (float("nan"), float("nan"))
    h = st.t.ppf(0.5 + conf / 2, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return (x.mean() - h, x.mean() + h)


def report(seeds, out_json=None):
    labels = {"val": per_mol("energy_largeval_base.json", "label"),
              "n19": per_mol("energy_n19_all.json", "label"),
              "small": per_mol("energy_small_baselines.json", "label")}
    names = {"val": [ln.strip() for ln in open(ROOT / "pretrain/rl/large_val.txt") if ln.strip()],
             "n19": [ln.strip() for ln in open(ROOT / "gauge_study/names_norb19.txt") if ln.strip()],
             "small": [ln.strip() for ln in open(ROOT / "pretrain/rl/small_val.txt") if ln.strip()]}
    logs = {s: read_log(s) for s in seeds}
    logs = {s: L for s, L in logs.items() if L and L["val"]}
    seeds = sorted(logs)
    summary = {"seeds": seeds, "per_seed": {}, "across": {}}

    # ---- val curves -----------------------------------------------------------------------------------------
    steps = sorted(set().union(*[set(L["val"]) for L in logs.values()]))
    lab_val = 100 * np.mean([labels["val"][n] for n in names["val"]])
    print(f"## Val curves (norb 17-18, 22 molecules, exact; label mean {lab_val:.2f})\n")
    print("| seed | " + " | ".join(str(s) for s in steps) + " | best step | done |")
    print("|:--|" + "--:|" * len(steps) + "--:|:--|")
    for s in seeds:
        v = logs[s]["val"]
        b = best_step(v)
        print(f"| {s} | " + " | ".join(f"{100 * v[t]['corr_mean']:.2f}" if t in v else "" for t in steps)
              + f" | {b} | {'yes' if logs[s]['done'] else 'running'} |")
    full = [s for s in seeds if all(t in logs[s]["val"] for t in steps)]
    rows = {"mean": [], "std": [], "min": [], "max": []}
    for t in steps:
        x = [100 * logs[s]["val"][t]["corr_mean"] for s in seeds if t in logs[s]["val"]]
        m, sd = mstd(x)
        rows["mean"].append(m)
        rows["std"].append(sd)
        rows["min"].append(min(x))
        rows["max"].append(max(x))
    n_at = [sum(t in logs[s]["val"] for s in seeds) for t in steps]
    print("| mean | " + " | ".join(f"{m:.2f}" for m in rows["mean"]) + " | | |")
    print("| std (ddof 1) | " + " | ".join(f"{x:.2f}" if np.isfinite(x) else "–" for x in rows["std"]) + " | | |")
    print("| range | " + " | ".join(f"{a:.1f}–{b:.1f}" for a, b in zip(rows["min"], rows["max"])) + " | | |")
    print("| seeds | " + " | ".join(str(k) for k in n_at) + " | | |")
    print()
    summary["curve"] = {"steps": steps, "label_mean": lab_val, **rows, "n_seeds": n_at,
                        "per_seed": {s: {t: 100 * logs[s]["val"][t]["corr_mean"] for t in logs[s]["val"]} for s in seeds},
                        "per_seed_median": {s: {t: 100 * logs[s]["val"][t]["corr_median"] for t in logs[s]["val"]}
                                            for s in seeds}}

    # ---- per seed: best / final on val (from the log), held-out norb 19 and small val --------------------------
    print("## Per seed: best (val-selected) and final policy\n")
    print("| seed | best step | val@best | gain | better | val@50 | gain | better | n19 best | gain | better | "
          "n19 final | gain | better | small best | small final | train 0→50 |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for s in seeds:
        L = logs[s]
        v = L["val"]
        b = best_step(v)
        rec = {"best_step": b}
        for tag, t in (("best", b), ("final", max(v) if L["done"] else None)):
            if t is None:
                continue
            pmv = v[t]["per_mol"]
            st_ = stats(pmv, labels["val"], names["val"])
            st_["mean"] = 100 * v[t]["corr_mean"]                       # exact mean (per_mol is rounded to 1e-4)
            st_["gain"] = st_["mean"] - lab_val
            st_["step"] = t
            rec[f"val_{tag}"] = st_
        for kind in ("n19", "small"):
            for which in ("best", "last"):
                pm = held_out(s, which, kind)
                if pm and all(n in pm for n in names[kind]):
                    rec[f"{kind}_{which}"] = stats(pm, labels[kind], names[kind])
                    if kind == "n19":
                        rec[f"{kind}_{which}_distinct"] = stats(pm, labels[kind], distinct(names[kind]))
        te = L["train_eval"]
        rec["train"] = {t: 100 * te[t]["corr_mean"] for t in te}
        summary["per_seed"][s] = rec

        def f(k, fld="mean", fmt="{:.1f}"):
            return fmt.format(rec[k][fld]) if k in rec else "–"

        def fb(k):
            return f"{rec[k]['better']}/{rec[k]['n']}" if k in rec else "–"
        tr = " → ".join(f"{x:.1f}" for _, x in sorted(rec["train"].items()))
        print(f"| {s} | {b} | {f('val_best', fmt='{:.2f}')} | {f('val_best', 'gain', '{:+.2f}')} | {fb('val_best')} | "
              f"{f('val_final', fmt='{:.2f}')} | {f('val_final', 'gain', '{:+.2f}')} | {fb('val_final')} | "
              f"{f('n19_best', fmt='{:.2f}')} | {f('n19_best', 'gain', '{:+.2f}')} | {fb('n19_best')} | "
              f"{f('n19_last', fmt='{:.2f}')} | {f('n19_last', 'gain', '{:+.2f}')} | {fb('n19_last')} | "
              f"{f('small_best', fmt='{:.2f}')} | {f('small_last', fmt='{:.2f}')} | {tr} |")
    print()

    # ---- across seeds ----------------------------------------------------------------------------------------
    print("## Across seeds (mean ± sample std over seeds; t 95% CI of the mean; per-molecule SE of one seed)\n")
    print("| quantity | seeds | values | mean ± std | 95% CI of mean | min–max | typical per-molecule SE |")
    print("|:--|--:|:--|--:|--:|--:|--:|")
    q = [("val@best (selected on this set)", "val_best", "mean"), ("val@best gain vs label", "val_best", "gain"),
         ("val@50", "val_final", "mean"), ("val@50 gain vs label", "val_final", "gain"),
         ("norb 19 best", "n19_best", "mean"), ("norb 19 best gain vs label", "n19_best", "gain"),
         ("norb 19 final", "n19_last", "mean"), ("norb 19 final gain vs label", "n19_last", "gain"),
         ("norb 19 best gain, 9 distinct", "n19_best_distinct", "gain"),
         ("norb 19 final gain, 9 distinct", "n19_last_distinct", "gain"),
         ("small val best", "small_best", "mean"), ("small val final", "small_last", "mean"),
         ("small val best gain vs label", "small_best", "gain"), ("small val final gain vs label", "small_last", "gain")]
    for title, k, fld in q:
        ss = [s for s in seeds if k in summary["per_seed"][s]]
        if not ss:
            continue
        x = [summary["per_seed"][s][k][fld] for s in ss]
        m, sd = mstd(x)
        lo, hi = tci(x)
        se = (np.mean([summary["per_seed"][s][k]["gain_se_mols"] for s in ss]) if fld == "gain" else float("nan"))
        summary["across"][f"{k}.{fld}"] = {"seeds": ss, "values": x, "mean": m, "std": sd, "ci95": [lo, hi],
                                           "per_mol_se": se}
        print(f"| {title} | {len(ss)} | {', '.join(f'{v:.2f}' for v in x)} | {m:.2f} ± {sd:.2f} | "
              f"{lo:.2f} – {hi:.2f} | {min(x):.2f}–{max(x):.2f} | {'' if not np.isfinite(se) else f'{se:.2f}'} |")
    print()
    # within-run wander of the val curve (eval-to-eval), for comparison with the between-seed spread
    wander = {}
    for s in seeds:
        v = logs[s]["val"]
        late = [100 * v[t]["corr_mean"] for t in sorted(v) if t >= 20]
        if len(late) >= 3:
            wander[s] = float(np.std(late, ddof=1))
    if wander:
        print("Within-run eval-to-eval std of val mean over steps 20-50: "
              + ", ".join(f"seed {s} {w:.2f}" for s, w in wander.items()) + "\n")
    summary["within_run_std_steps20_50"] = wander
    summary["extras"] = extras(summary, logs, labels, names, seeds, lab_val)
    if out_json:
        Path(out_json).parent.mkdir(parents=True, exist_ok=True)
        json.dump(summary, open(out_json, "w"), indent=1, default=float)
        print(f"-> {out_json}")
    return summary


def extras(summary, logs, labels, names, seeds, lab_val):
    """Added on resume (5 Oct): the shared start policy (rl4s = grpo_slot_T4p best, step 0 of every seed) on the
    held-out sets, each seed's mean val over steps 20-50 (less selection-biased than the best step), and per-molecule
    agreement of the seeds on the val set."""
    ex = {}
    init = {"n19": per_mol("energy_n19_all.json", "rl4s"), "small": per_mol("energy_smallval_rl4L.json", "rl4s")}
    print("## Extras\n")
    print("| set | label | start policy rl4s (step 0) | seed: best − rl4s | seed: final − rl4s |")
    print("|:--|--:|--:|:--|:--|")
    for kind, title in (("n19", "norb 19 (12)"), ("small", "norb 15-16 small val (9)")):
        if not (init[kind] and all(n in init[kind] for n in names[kind])):
            continue
        st0 = stats(init[kind], labels[kind], names[kind])
        ex[f"{kind}_init"] = st0
        lab = st0["mean"] - st0["gain"]
        d = {w: {s: summary["per_seed"][s][f"{kind}_{w}"]["mean"] - st0["mean"] for s in seeds
                 if f"{kind}_{w}" in summary["per_seed"][s]} for w in ("best", "last")}
        ex[f"{kind}_minus_init"] = d
        print(f"| {title} | {lab:.2f} | {st0['mean']:.2f} ({st0['gain']:+.2f} vs label) | "
              + ", ".join(f"s{s} {v:+.2f}" for s, v in d["best"].items()) + " | "
              + ", ".join(f"s{s} {v:+.2f}" for s, v in d["last"].items()) + " |")
    print()
    late = {s: float(np.mean([100 * logs[s]["val"][t]["corr_mean"] for t in sorted(logs[s]["val"]) if t >= 20]))
            for s in seeds}
    ex["val_mean_steps20_50"] = late
    m, sd = mstd(list(late.values()))
    print("Mean val over steps 20-50 per seed (gain vs label): "
          + ", ".join(f"seed {s} {v:.2f} ({v - lab_val:+.2f})" for s, v in late.items())
          + f"; across seeds {m:.2f} ± {sd:.2f}\n")
    for tag in ("best", "final"):
        recs = [summary["per_seed"][s].get(f"val_{tag}") for s in seeds]
        if any(r is None for r in recs):
            continue
        P = np.array([[100 * logs[s]["val"][r["step"]]["per_mol"][n] for n in names["val"]] for s, r in zip(seeds, recs)])
        lab = 100 * np.array([labels["val"][n] for n in names["val"]])
        sd_m = P.std(0, ddof=1)
        allw, alll = int((P > lab).all(0).sum()), int((P < lab).all(0).sum())
        ex[f"val_{tag}_per_mol"] = {"all_seeds_beat_label": allw, "all_seeds_below_label": alll,
                                    "per_mol_sd_over_seeds_mean": float(sd_m.mean()),
                                    "per_mol_sd_over_seeds_max": float(sd_m.max())}
        print(f"Val {tag} step, per molecule over {len(seeds)} seeds: all seeds beat the label on {allw}/22, all below on "
              f"{alll}/22; per-molecule SD over seeds mean {sd_m.mean():.2f}, max {sd_m.max():.2f} points\n")
    return ex


SEED_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]     # dataviz reference palette, slots 1-4 (validated)
SEED_MARKERS = ["o", "s", "^", "D"]
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"


def figure(summary, out_png):
    """(a) val curve per seed + across-seed mean and std band; (b) gain over the labels per seed and set."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#c3c2b7",
                         "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED})
    cur = summary["curve"]
    seeds = summary["seeds"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw={"width_ratios": [1.35, 1]}, facecolor=SURF)
    for ax in (a, b):
        ax.set_facecolor(SURF)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    steps = np.array(cur["steps"])
    m, sd = np.array(cur["mean"]), np.array(cur["std"])
    ok = np.isfinite(sd) & (np.array(cur["n_seeds"]) == len(seeds))      # only steps every seed has reached
    if ok.sum() > 1:
        a.fill_between(steps[ok], (m - sd)[ok], (m + sd)[ok], color="#c3c2b7", alpha=0.45, lw=0,
                       label="mean ± std over seeds")
    for i, s in enumerate(seeds):
        ps = cur["per_seed"][s]
        t = sorted(ps)
        a.plot(t, [ps[k] for k in t], color=SEED_COLORS[i % 4], lw=2, marker=SEED_MARKERS[i % 4], ms=6.5,
               mec=SURF, mew=1.2, label=f"seed {s}" + (" (documented run)" if s == 0 else ""))
    if ok.sum() > 1 and len(seeds) > 1:
        a.plot(steps[ok], m[ok], color=INK, lw=2.4, label="mean over seeds")
    a.axhline(cur["label_mean"], color=MUTED, lw=1.5, ls=(0, (4, 3)))
    a.text(steps[0], cur["label_mean"] + 0.2, f"optimize=True labels {cur['label_mean']:.1f}", color=INK2, ha="left",
           va="bottom", fontsize=9)
    a.set_xlabel("GRPO step")
    a.set_ylabel("val mean % of CCSD corr. energy")
    a.set_title("(a) norb 17–18 val (22 molecules, exact)", loc="left", color=INK, fontsize=11)
    a.set_xticks(steps)
    a.legend(frameon=False, fontsize=8.5, loc="lower right", ncol=2, bbox_to_anchor=(1.0, 0.07))
    # (b) gains over the labels
    cats = [("val_best", "17–18 val\nbest step*"), ("val_final", "17–18 val\nstep 50"), ("n19_best", "norb 19\nbest step"),
            ("n19_last", "norb 19\nstep 50"), ("small_best", "15–16 val\nbest step"), ("small_last", "15–16 val\nstep 50")]
    for j, (k, _) in enumerate(cats):
        xs = []
        for i, s in enumerate(seeds):
            r = summary["per_seed"][s].get(k)
            if r is None:
                continue
            x = j + (i - (len(seeds) - 1) / 2) * 0.13
            b.plot([x], [r["gain"]], marker=SEED_MARKERS[i % 4], color=SEED_COLORS[i % 4], ms=7.5, mec=SURF, mew=1.2,
                   ls="none", label=f"seed {s}" if j == 0 else None)
            xs.append(r["gain"])
        if len(xs) > 1:
            b.plot([j - 0.3, j + 0.3], [np.mean(xs)] * 2, color=INK, lw=2.2, label="mean" if j == 0 else None)
    b.axhline(0, color="#c3c2b7", lw=1.2)
    b.set_xticks(range(len(cats)))
    b.set_xticklabels([c[1] for c in cats], fontsize=8.5)
    b.set_ylabel("policy − label (% points)")
    b.set_title("(b) gain over the optimize=True labels", loc="left", color=INK, fontsize=11)
    b.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    fig.text(0.995, 0.01, "* selected on the same val set (optimistic)", ha="right", va="bottom", fontsize=8, color=MUTED)
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, facecolor=SURF)
    print(f"-> {out_png}")


SETS = {"n19": ("gauge_study/names_norb19.txt", "n19"), "small": ("pretrain/rl/small_val.txt", "smallval")}


def eval_cmd(seed, kind, which, queue_root, threads):
    """Energies of seed's best / last policy on a held-out set through the shared exact queue."""
    import os
    import subprocess
    import sys

    import torch
    names_file, short = SETS[kind]
    d = run_dir(seed)
    steps = {w: torch.load(d / f"policy_{w}.pt", map_location="cpu", weights_only=False)["step"] for w in which}
    todo, alias = [], {}
    for w in which:
        twin = [v for v in todo if steps[v] == steps[w]]
        if twin:
            alias[w] = twin[0]
        else:
            todo.append(w)
    env = dict(os.environ, **{v: str(threads) for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                         "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS")})
    dump = f"runs_ot/energy_tasks/task7_{short}_seed{seed}.pkl"
    out = FRES / f"energy_{short}_seed{seed}.json"
    if out.exists():
        raise SystemExit(f"{out} exists")
    cmd = [sys.executable, "-m", "pretrain.rl.policy_dump", "--names-file", names_file, "--device", "cpu", "--out", dump]
    for w in todo:
        cmd += ["--policy", str((d / f"policy_{w}.pt").relative_to(ROOT)), "--tag", f"s{seed}{w}"]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT, env=env)
    cmd = [sys.executable, "-m", "pretrain.rl.queue_eval", "--dump", dump, "--queue-root", queue_root,
           "--out", str(out.relative_to(ROOT))]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT, env=env)
    res = json.load(open(out))
    res["policy_steps"] = {f"s{seed}{w}": steps[w] for w in which}
    for w, v in alias.items():                       # same checkpoint: store the result under both tags
        for n in res["per_molecule"]:
            res["per_molecule"][n][f"s{seed}{w}"] = dict(res["per_molecule"][n][f"s{seed}{v}"], same_as=f"s{seed}{v}")
    json.dump(res, open(out, "w"), indent=1)
    print(f"steps {steps}; aliases {alias} -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("report")
    r.add_argument("--seeds", type=int, nargs="*", default=[0, 1, 2, 3])
    r.add_argument("--json", default=str(FRES / "summary.json"))
    r.add_argument("--png", default=str(ROOT / "docs" / "followups" / "task7_rl_seeds.png"))
    e = sub.add_parser("eval")
    e.add_argument("--seed", type=int, required=True)
    e.add_argument("--kind", choices=sorted(SETS), required=True)
    e.add_argument("--which", nargs="+", default=["best", "last"], choices=["best", "last"])
    e.add_argument("--queue-root", default="rl_queue/exact_shared")
    e.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    if a.cmd == "report":
        summ = report(a.seeds, a.json)
        if a.png:
            figure(summ, a.png)
    elif a.cmd == "eval":
        eval_cmd(a.seed, a.kind, a.which, a.queue_root, a.threads)


if __name__ == "__main__":
    main()
