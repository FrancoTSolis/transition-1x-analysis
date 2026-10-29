"""Figures of the UCJ-init study (docs/ucj_init/).

  fig_bch2_runaway.png   BCH2 energy per cycle of Avijit's optimizer (unmodified): small systems (his default seed)
                         and a TS molecule from five seeds; CCSD correlation energy for scale.
  fig_exact_opt.png      exact-energy L-BFGS from each start: best-so-far % of the CCSD correlation energy vs
                         gradient evaluations, one panel per TS molecule.

    python -m ucj_init.make_figs
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from ucj_init.report_exact import load

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "pretrain" / "opt_true" / "results" / "ucj_init"
OUT = ROOT / "docs" / "ucj_init"

# reference palette slots 1-3 (validated all-pairs), slot 7 with a dash pattern, gray for naive starts
BLUE, ORANGE, AQUA, VIOLET, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#8c8b86"
INK, MUTED = "#0b0b0b", "#52514e"
STYLE = {"rl4n29f": (BLUE, "-", "network + RL (rl4n29f)"),
         "label": (VIOLET, (0, (5, 2)), "our optimize=True label"),
         "truncated": (ORANGE, "-", "Lin: truncated t2 (optimize=False)"),
         "cdf_lin": (AQUA, (0, (2, 1.5)), "Lin: compressed DF (optimize=True)"),
         "random0": (GRAY, "-", "random")}


def _axes_style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(alpha=0.25, lw=0.6)


def fig_bch2():
    small = json.load(open(RES / "bch2_reproduce.json"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    ax = axes[0]
    cols = [BLUE, ORANGE, AQUA, VIOLET, GRAY, "#b0afa9", "#5f5e5a"]
    for c, r in zip(cols, small):
        y = [-h["e_corr"] for h in r["run_history"]]
        ax.plot(np.arange(1, len(y) + 1), y, color=c, lw=1.6, label=r["name"])
        ax.scatter([len(y)], [-r["ccsd_corr"]], color=c, s=14, marker="_")
    ax.set_yscale("log")
    ax.set_xlabel("optimizer cycle", color=INK, fontsize=9)
    ax.set_ylabel("$-(E_\\mathrm{BCH2} - E_\\mathrm{HF})$ (Eh)", color=INK, fontsize=9)
    ax.set_title("Avijit's default seed, small systems", fontsize=9, color=INK)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    _axes_style(ax)
    ax = axes[1]
    seeds = [("default", GRAY, "his default seed"), ("random", "#b0afa9", "random"),
             ("truncated", ORANGE, "Lin truncated"), ("label", VIOLET, "our label"), ("rl4n29f", BLUE, "network")]
    ccsd = None
    for seed, c, lab in seeds:
        f = RES / "bch2_ts" / f"C2H3N_rxn2858_TS__{seed}.json"
        if not f.exists():
            continue
        r = json.load(open(f))
        ccsd = -r["ccsd_corr"]
        y = [max(1e-6, -h["e_corr"]) for h in r["history"]]
        ax.plot(np.arange(1, len(y) + 1), y, color=c, lw=1.6, label=lab)
    if ccsd:
        ax.axhline(ccsd, color=MUTED, lw=0.8, ls=":")
        ax.text(1, ccsd * 1.25, "CCSD correlation energy", fontsize=7, color=MUTED)
    ax.set_yscale("log")
    ax.set_xlabel("optimizer cycle", color=INK, fontsize=9)
    ax.set_title("C2H3N_rxn2858_TS (norb 15), five seeds", fontsize=9, color=INK)
    ax.legend(fontsize=7, frameon=False)
    _axes_style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig_bch2_runaway.png", dpi=140)


def fig_exact(method="lbfgs"):
    runs = load(method)
    names = sorted({n for n, _ in runs}, key=lambda n: (n.split("_")[0], n))
    names = [n for n in names if (n, "rl4n29f") in runs and (n, "truncated") in runs
             and (method != "lbfgs" or n.endswith("_TS"))]
    ncol = min(3, len(names))
    nrow = (len(names) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), squeeze=False)
    for ax, name in zip(axes.flat, names):
        for init, (c, ls, lab) in STYLE.items():
            r = runs.get((name, init))
            if r is None:
                continue
            b = r["best"]
            x = np.arange(len(b))
            ax.plot(x, np.where(np.isfinite(b), b, np.nan), color=c, ls=ls, lw=1.7, label=lab)
        ax.set_xscale("symlog", linthresh=10)
        ax.set_ylim(0, 95)
        ax.set_title(name, fontsize=9, color=INK)
        ax.set_xlabel("gradient evaluations", color=INK, fontsize=8)
        ax.set_ylabel("% of CCSD corr. (best so far)", color=INK, fontsize=8)
        _axes_style(ax)
    # one legend for every start that appears in any panel (the first panel may lack some)
    used = [i for i in STYLE if any((n, i) in runs for n in names)]
    handles = [Line2D([], [], color=STYLE[i][0], ls=STYLE[i][1], lw=1.7, label=STYLE[i][2]) for i in used]
    spare = list(axes.flat)[len(names):]
    for ax in spare:
        ax.axis("off")
    if spare:
        spare[0].legend(handles=handles, fontsize=9, frameon=False, loc="center")
        fig.tight_layout()
    else:
        fig.legend(handles=handles, fontsize=8, frameon=False, loc="upper center", ncol=len(handles))
        fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / f"fig_exact_opt_{method}.png", dpi=140)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_bch2()
    fig_exact("lbfgs")
    if any(Path(ROOT / "rl_runs" / "ucj_init" / "exact_opt").glob("*__trn.jsonl")):
        fig_exact("trn")
    print("figures in", OUT)
