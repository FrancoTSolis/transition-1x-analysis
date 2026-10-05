#!/usr/bin/env python3
"""Figures for docs/followups/task2_direct_opt.md (static PNG; the tables carry every number).

  fig_task2_lucj.png : objective A, best-so-far % corr vs evaluations, one panel per molecule; color = optimizer,
                       line style = start; dotted = the one-call network value
  fig_task2_long.png : the two 5000-evaluation SPSA runs next to the 500-evaluation ones
  fig_task2_qsci.png : objective B, best-so-far QSCI % corr (top) and the subspace size of each evaluation
                       (bottom, strings per spin), one column per molecule
Palette: categorical slots 1-3 of the dataviz reference palette (blue, orange, aqua), fixed order.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.followups import task2_common as C  # noqa: E402

S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
OUT = C.ROOT / "docs" / "followups"


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def pct_trace(r):
    return [C.corr_pct(v, r["e_hf"], r["e_ccsd"]) for v in r["best_trace"]]


def fig_lucj():
    """One panel per molecule (3 columns): best-so-far exact LUCJ energy of the 500-evaluation runs.  Color = optimizer
    (NOMAD blue, SPSA orange), line style = start (solid: label, dashed: network + RL rl4n29f, SPSA with the label-start
    gain); dotted gray = the network's one-call value."""
    res = {}
    for f in sorted((C.RESULTS / "lucj").glob("*.json")):
        r = json.load(open(f))
        res[(r["name"], r["start"], r["optimizer"], r["args"].get("tag") or "")] = r
    names = sorted({k[0] for k in res if k[3] in ("", "aLabel")}, key=lambda n: (C.load_ham(n)["norb"], n))
    names = [n for n in names if (n, "label", "spsa", "") in res or (n, "label", "nomad", "") in res]
    if not names:
        return
    ncol = 3
    nrow = (len(names) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.7 * nrow + 0.6), sharex=True, squeeze=False)
    fig.patch.set_facecolor(SURF)
    series = ((("label", "nomad", ""), S1, "-", "NOMAD from label"),
              (("label", "spsa", ""), S2, "-", "SPSA from label"),
              (("rl4n29f", "nomad", ""), S1, (0, (5, 2.5)), "NOMAD from network + RL"),
              (("rl4n29f", "spsa", "aLabel"), S2, (0, (5, 2.5)), "SPSA from network + RL (label-start gain)"))
    for j, n in enumerate(names):
        ax = axes[j // ncol, j % ncol]
        style(ax)
        rl = next((res[(n, "rl4n29f", o, t)] for o, t in (("nomad", ""), ("spsa", ""), ("spsa", "aLabel"))
                   if (n, "rl4n29f", o, t) in res), None)
        for (st, o, tag), col, ls, lab in series:
            r = res.get((n, st, o, tag))
            if r:
                y = pct_trace(r)
                ax.plot(np.arange(1, len(y) + 1), y, color=col, lw=2, ls=ls, solid_capstyle="round", label=lab)
        if rl is not None:
            ax.axhline(rl["corr_var0"], color=INK2, lw=1, ls=(0, (1, 2)), label="network + RL, one call")
        ax.set_title(f"{n} (norb {C.load_ham(n)['norb']})", fontsize=8.5, color=INK, loc="left")
        if j % ncol == 0:
            ax.set_ylabel("% of CCSD corr. (best so far)", fontsize=8, color=INK2)
        if j // ncol == nrow - 1:
            ax.set_xlabel("exact LUCJ energy evaluations", fontsize=8, color=INK2)
    for j in range(len(names), nrow * ncol):
        axes[j // ncol, j % ncol].set_visible(False)
    h, l = [], []
    for ax in axes.ravel():
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh)
                l.append(ll)
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8, labelcolor=INK,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Objective A: exact LUCJ energy, 500 evaluations per molecule and start", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(OUT / "fig_task2_lucj.png", dpi=150, bbox_inches="tight", facecolor=SURF)
    plt.close(fig)
    print("->", OUT / "fig_task2_lucj.png")


def fig_long():
    """The 5000-evaluation SPSA runs next to the 500-evaluation ones (same molecule), with the one-call value.  Same
    encoding as fig_task2_lucj (SPSA orange; solid label start, dashed network start); thin = 500-evaluation runs."""
    runs = {}
    for f in sorted((C.RESULTS / "lucj").glob("*__spsa*.json")):
        r = json.load(open(f))
        tag = r["args"].get("tag") or ""
        runs[(r["name"], r["start"], tag)] = r
    longs = [k for k in runs if k[2].startswith("long")]
    if not longs:
        return
    n = longs[0][0]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    fig.patch.set_facecolor(SURF)
    style(ax)
    dash = (0, (5, 2.5))
    for (st, tag), ls, lw, lab in ((("label", "long5000"), "-", 2.2, "SPSA from label, 5000 evaluations"),
                                   (("rl4n29f", "long5000aLabel"), dash, 2.2,
                                    "SPSA from network + RL, 5000 evaluations"),
                                   (("label", ""), "-", 1.0, "SPSA from label, 500-evaluation run"),
                                   (("rl4n29f", "aLabel"), dash, 1.0, "SPSA from network + RL, 500-evaluation run")):
        r = runs.get((n, st, tag))
        if r:
            y = pct_trace(r)
            ax.plot(np.arange(1, len(y) + 1), y, color=S2, lw=lw, ls=ls, label=lab)
    r1 = runs.get((n, "rl4n29f", "aLabel")) or runs.get((n, "rl4n29f", "long5000aLabel"))
    if r1:
        ax.axhline(r1["corr_var0"], color=INK2, lw=1, ls=(0, (1, 2)), label="network + RL, one call")
    ax.set_xlabel("exact LUCJ energy evaluations", fontsize=8, color=INK2)
    ax.set_ylabel("% of CCSD corr. (best so far)", fontsize=8, color=INK2)
    ax.set_title(f"{n}: long SPSA runs (gain schedule set by the budget)", fontsize=9, color=INK, loc="left")
    ax.legend(loc="lower right", frameon=False, fontsize=7.5, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(OUT / "fig_task2_long.png", dpi=150, bbox_inches="tight", facecolor=SURF)
    plt.close(fig)
    print("->", OUT / "fig_task2_long.png")


def fig_qsci():
    runs = []
    for f in sorted((C.LOGS / "qsci").glob("*__nomad.jsonl")) + sorted((C.LOGS / "qsci").glob("*__spsa.jsonl")):
        name, start, opt = f.stem.split("__")[:3]
        h = [json.loads(ln) for ln in open(f)]
        if len(h) < 5:
            continue
        runs.append((name, start, opt, h))
    names = sorted({r[0] for r in runs})
    if not names:
        return
    fci_p = C.ROOT / "pretrain/opt_true/results/followups/task1_sqd_vs_lin/fci_refs.json"
    fci = json.load(open(fci_p)) if fci_p.exists() else {}
    fig, axes = plt.subplots(2, len(names), figsize=(4.2 * len(names), 5.6), sharex="col", squeeze=False)
    fig.patch.set_facecolor(SURF)
    cols = {("label", "nomad"): (S1, "NOMAD from label"), ("label", "spsa"): (S2, "SPSA from label"),
            ("rl4n29f", "nomad"): (S3, "NOMAD from rl4n29f")}
    for j, n in enumerate(names):
        ham = C.load_ham(n)
        ax, bx = axes[0, j], axes[1, j]
        style(ax)
        style(bx)
        for name, start, opt, h in runs:
            if name != n or (start, opt) not in cols:
                continue
            col, lab = cols[(start, opt)]
            x = np.array([r["n"] for r in h])
            best = np.array([C.corr_pct(r["best"], ham["e_hf"], ham["e_ccsd"]) for r in h])
            ax.plot(x, best, color=col, lw=2, label=lab)
            bx.plot(x, [r["dim"][0] for r in h], color=col, lw=1, alpha=0.9, label=lab)
        if n in fci:
            yf = C.corr_pct(fci[n]["e_fci"], ham["e_hf"], ham["e_ccsd"])
            ax.axhline(yf, color=INK2, lw=1, ls=(0, (4, 3)))
            ax.text(ax.get_xlim()[1], yf, "FCI", color=INK2, fontsize=7, ha="right", va="bottom")
        ax.set_title(n, fontsize=8.5, color=INK, loc="left")
        ax.set_ylabel("QSCI of 10^4 samples, % CCSD corr.\n(best so far)", fontsize=8, color=INK2)
        bx.set_ylabel("strings per spin in the subspace", fontsize=8, color=INK2)
        bx.set_xlabel("QSCI evaluations", fontsize=8, color=INK2)
    h, l = [], []
    for ax in axes.ravel():
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh)
                l.append(ll)
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8, labelcolor=INK,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Objective B: Lin et al.'s QSCI objective (exact samples)", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(OUT / "fig_task2_qsci.png", dpi=150, bbox_inches="tight", facecolor=SURF)
    print("->", OUT / "fig_task2_qsci.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_lucj()
    fig_long()
    fig_qsci()
