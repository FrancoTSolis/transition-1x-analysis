"""Figure for slides/brandon_2026-10-29.md: exact-energy optimization of one TS molecule from each start,
L-BFGS and Avijit's trust-region Newton logic side by side.

    python slides/make_brandon_figs.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from ucj_init.make_figs import STYLE, _axes_style, INK
from ucj_init.report_exact import load

NAME = "C2H4O_rxn0724_TS"
OUT = Path(__file__).resolve().parent / "fig_brandon_convergence.png"


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
    for ax, (method, title) in zip(axes, [("lbfgs", "L-BFGS"), ("trn", "Avijit's trust-region Newton logic")]):
        runs = load(method)
        for init, (c, ls, lab) in STYLE.items():
            r = runs.get((NAME, init))
            if r is None:
                continue
            b = r["best"]
            ax.plot(np.arange(len(b)), np.where(np.isfinite(b), b, np.nan), color=c, ls=ls, lw=2.0)
        ax.set_xscale("symlog", linthresh=10)
        ax.set_ylim(0, 95)
        ax.set_title(f"{title}", fontsize=11, color=INK)
        ax.set_xlabel("gradient evaluations", color=INK, fontsize=10)
        _axes_style(ax)
    axes[0].set_ylabel("% of CCSD correlation energy", color=INK, fontsize=10)
    handles = [Line2D([], [], color=c, ls=ls, lw=2.0, label=lab) for c, ls, lab in STYLE.values()]
    fig.legend(handles=handles, fontsize=9, frameon=False, loc="upper center", ncol=len(handles))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
