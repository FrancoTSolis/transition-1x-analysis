"""Tables and figures of the exact-energy optimization runs (ucj_init.opt_exact).

Per (molecule, init, method): best-so-far % of the CCSD correlation energy after B gradient evaluations,
evaluations needed to reach thresholds, final value, and the evaluations a start needs to reach the network's
one-call value.

    python -m ucj_init.report_exact [--method lbfgs] [--fig docs/ucj_init/fig_exact_opt.png]
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / "rl_runs" / "ucj_init" / "exact_opt"
RES = ROOT / "pretrain" / "opt_true" / "results" / "ucj_init" / "exact_opt"
ORDER = ["zero", "random0", "random1", "truncated", "cdf_lin", "label", "pre4", "rl4L", "rl4n29f"]
BUDGETS = [0, 10, 25, 50, 100, 200, 300, 400]
THRESH = [60, 70, 75, 80, 85]


def load(method):
    runs = {}
    for f in sorted(LOGS.glob(f"*__{method}.jsonl")):
        name, init, _ = f.stem.split("__")
        recs = [json.loads(l) for l in open(f)]
        if not recs:
            continue
        # best so far as a function of the evaluation count (record "i" counts every gradient evaluation,
        # including the Hessian-vector products of the trust-region method, which are not logged)
        n_ev = recs[-1]["i"]
        pct = np.full(n_ev, -np.inf)
        for r in recs:
            pct[r["i"] - 1] = max(pct[r["i"] - 1], r["pct"])
        res = RES / f"{f.stem}.json"
        runs[(name, init)] = dict(pct=pct, best=np.maximum.accumulate(pct), gnorm=np.array([r["gnorm"] for r in recs]),
                                  t=np.array([r["t"] for r in recs]), done=res.exists(),
                                  meta=json.load(open(res)) if res.exists() else None)
    return runs


def at(best, b):
    return best[min(b, len(best) - 1)] if b < len(best) or len(best) else np.nan


def evals_to(best, thr):
    idx = np.nonzero(best >= thr)[0]
    return int(idx[0]) if len(idx) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", default="lbfgs")
    ap.add_argument("--fig", default=None)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    runs = load(args.method)
    names = sorted({n for n, _ in runs})
    inits = [i for i in ORDER if any((n, i) in runs for n in names)]
    summary = {}
    print(f"method {args.method}; % of CCSD correlation energy, best so far after B evaluations "
          f"(eval 0 = the start); * = run still going\n")
    for name in names:
        print(f"## {name}")
        print(f"{'init':10s} " + " ".join(f"B={b:<4d}" for b in BUDGETS) + "  n_eval  final|g|  " +
              " ".join(f">={t}%" for t in THRESH))
        for init in inits:
            r = runs.get((name, init))
            if r is None:
                continue
            best = r["best"]
            row = [f"{best[b]:6.1f}" if b < len(best) else "     -" for b in BUDGETS]
            ev = [evals_to(best, t) for t in THRESH]
            print(f"{init:10s} " + " ".join(row) + f"  {len(best):5d}{'*' if not r['done'] else ' '}  {r['gnorm'][-1]:.1e}  "
                  + " ".join(f"{('-' if e is None else e):>5}" for e in ev))
            summary[f"{name}|{init}"] = dict(budgets={b: float(best[b]) for b in BUDGETS if b < len(best)},
                                             evals_to={t: e for t, e in zip(THRESH, ev)}, n=len(best),
                                             final=float(best[-1]), done=r["done"])
        # evaluations each start needs to reach the network's one-call value
        if (name, "rl4n29f") in runs:
            target = runs[(name, "rl4n29f")]["best"][0]
            eq = {i: evals_to(runs[(name, i)]["best"], target) for i in inits if (name, i) in runs}
            print(f"evals to reach the one-call value of rl4n29f ({target:.1f}%): " +
                  ", ".join(f"{i} {('>' + str(len(runs[(name, i)]['best']))) if e is None else e}" for i, e in eq.items()))
        print()
    # aggregate over molecules that have every init
    print("## mean over molecules (those with the init)")
    print(f"{'init':10s} " + " ".join(f"B={b:<4d}" for b in BUDGETS) + "   n_mol")
    for init in inits:
        rows = [runs[(n, init)]["best"] for n in names if (n, init) in runs]
        vals = []
        for b in BUDGETS:
            v = [r[b] for r in rows if b < len(r)]
            vals.append(f"{np.mean(v):6.1f}" if len(v) == len(rows) and v else "     -")
        print(f"{init:10s} " + " ".join(vals) + f"   {len(rows)}")
    if args.json:
        json.dump(summary, open(args.json, "w"), indent=1)
    if args.fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cols = {"zero": "#999999", "random0": "#bbbbbb", "random1": "#cccccc", "truncated": "#e08214",
                "cdf_lin": "#b35806", "label": "#8073ac", "pre4": "#4393c3", "rl4L": "#2166ac", "rl4n29f": "#053061"}
        nn = len(names)
        ncol = min(3, nn)
        nrow = (nn + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.6 * nrow), squeeze=False)
        for ax, name in zip(axes.flat, names):
            for init in inits:
                r = runs.get((name, init))
                if r is None:
                    continue
                ax.plot(np.arange(len(r["best"])), r["best"], color=cols.get(init, "k"), lw=1.6, label=init)
            ax.set_title(name, fontsize=9)
            ax.set_xlabel("gradient evaluations")
            ax.set_ylabel("% CCSD corr. (best so far)")
            ax.set_ylim(-20, 100)
            ax.set_xscale("symlog", linthresh=10)
            ax.grid(alpha=0.3)
        for ax in list(axes.flat)[nn:]:
            ax.axis("off")
        axes.flat[0].legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        Path(args.fig).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.fig, dpi=130)
        print("figure:", args.fig)


if __name__ == "__main__":
    main()


def paired_table(method="lbfgs", ref="rl4n29f", budgets=(0, 10, 25, 50, 100, 150), names_filter=None):
    """Markdown table: each start against `ref` on the molecules where both ran (budgets both runs reached)."""
    runs = load(method)
    names = sorted({n for n, _ in runs if names_filter is None or names_filter(n)})
    lines = ["| start | molecules | " + " | ".join(f"% at {b}" for b in budgets) +
             " | evals to 70 % (median) | evals to 80 % (median) | network ahead at 100 |",
             "|:--|--:|" + "--:|" * len(budgets) + "--:|--:|--:|"]
    for init in ORDER:
        mols = [n for n in names if (n, init) in runs and (n, ref) in runs]
        if not mols:
            continue
        cells = []
        for b in budgets:
            v = [runs[(n, init)]["best"][b] for n in mols if b < len(runs[(n, init)]["best"])]
            cells.append(f"{np.mean(v):.1f}" if len(v) == len(mols) else "–")

        def med(thr):
            e = [evals_to(runs[(n, init)]["best"], thr) for n in mols]
            e = [x if x is not None else np.inf for x in e]
            m = float(np.median(e))
            return "never" if not np.isfinite(m) else f"{m:.0f}"
        ahead = [runs[(n, ref)]["best"][100] > runs[(n, init)]["best"][100] for n in mols
                 if len(runs[(n, ref)]["best"]) > 100 and len(runs[(n, init)]["best"]) > 100]
        lines.append(f"| {init} | {len(mols)} | " + " | ".join(cells) + f" | {med(70)} | {med(80)} | "
                     f"{'–' if init == ref else f'{sum(ahead)}/{len(ahead)}'} |")
    return "\n".join(lines)


def final_table(method="lbfgsJ", inits=("truncated", "cdf_lin", "label", "pre4", "rl4n29f")):
    """Markdown table: final best % per molecule (rows) and start (columns)."""
    runs = load(method)
    names = sorted({n for n, _ in runs})
    lines = ["| molecule | " + " | ".join(inits) + " |", "|:--|" + "--:|" * len(inits)]
    for n in names:
        lines.append(f"| {n} | " + " | ".join(f"{runs[(n, i)]['best'][-1]:.1f}" if (n, i) in runs else ""
                                            for i in inits) + " |")
    return "\n".join(lines)
