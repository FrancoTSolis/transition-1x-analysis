#!/usr/bin/env python3
"""Figure for docs/followups/task1_sqd_vs_lin.md: QSCI error against the LUCJ variational error and against the
subspace size (strings per spin), norb 15-16 small-val set, FCI reference.

  python3 -m pretrain.followups.task1_sqd_figure [--proto lin] [--out docs/followups/task1_sqd_vs_lin.png]

Gray dots: one (molecule, candidate) each (QSCI = mean over the independent 10^5-sample draws).  Blue markers:
candidate medians over the molecules, labelled directly with leader lines (one series, so no legend box).
"""
from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator  # noqa: E402

from pretrain.followups.task1_sqd_report import CANDS, LAB, ROOT, SETS, build, has, load_all, q_seeds  # noqa: E402

SURF, INK, INK2, MUTED, GRID, AXIS, BLUE = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#2a78d6"
# label anchor of each candidate's median marker, in axes-fraction coordinates (panel A, panel B)
POS = {"truncated": ((0.47, 0.90), (0.20, 0.80)),
       "cdf_lin": ((0.66, 0.42), (0.80, 0.45)),
       "label": ((0.30, 0.72), (0.42, 0.72)),
       "pre4": ((0.42, 0.57), (0.14, 0.60)),
       "rl4L": ((0.30, 0.08), (0.62, 0.07)),
       "rl4n29f": ((0.03, 0.60), (0.22, 0.20))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/followups/task1_sqd_vs_lin.png")
    ap.add_argument("--proto", default="lin")
    a = ap.parse_args()
    names = [ln.strip() for ln in open(ROOT / SETS["smallval"]) if ln.strip()]
    st, diag, fci, cct = load_all()
    rows = build(names, st, diag, fci, cct)
    pts = {c: [] for c in CANDS}
    for n in names:
        for c in CANDS:
            r = rows.get((n, c))
            if not has(r, a.proto) or r["refs"]["fci"] is None:
                continue
            R = r["refs"]["fci"]
            q = q_seeds(r, a.proto).mean()
            dim = np.mean([np.mean(r["seeds"][s][a.proto]["dim"]) for s in r["seeds"] if a.proto in r["seeds"][s]])
            pts[c].append((1e3 * (r["e_lucj"] - R), 1e3 * (q - R), dim))
    nmol = max(len(v) for v in pts.values())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": AXIS,
                         "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.titlecolor": INK, "axes.titlesize": 10})
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True, facecolor=SURF)
    yt = [2, 5, 10, 20, 50, 100]
    for k, (ax, xi, xl, xt) in enumerate([
            (axs[0], 0, "LUCJ variational energy − E_FCI (mHa)", [50, 100, 200, 500, 1000, 2000]),
            (axs[1], 2, "strings per spin in the QSCI subspace", [50, 100, 200, 500, 1000])]):
        ax.set_facecolor(SURF)
        ax.set_xscale("log")
        ax.set_yscale("log")
        for c in CANDS:
            if not pts[c]:
                continue
            p = np.array(pts[c])
            ax.scatter(p[:, xi], p[:, 1], s=14, color=MUTED, alpha=0.5, linewidths=0, zorder=2)
            mx, my = np.median(p[:, xi]), np.median(p[:, 1])
            ax.scatter([mx], [my], s=64, color=BLUE, edgecolors=SURF, linewidths=2, zorder=4)
            ax.annotate(LAB[c], xy=(mx, my), xytext=POS[c][k], textcoords="axes fraction", color=INK,
                        fontsize=8.5, zorder=5, va="center",
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=2, shrinkB=5))
        ax.xaxis.set_major_locator(FixedLocator(xt))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="major", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
        ax.set_xlabel(xl)
    axs[0].yaxis.set_major_locator(FixedLocator(yt))
    axs[0].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    axs[0].yaxis.set_minor_locator(NullLocator())
    axs[0].yaxis.set_minor_formatter(NullFormatter())
    axs[0].set_ylabel("QSCI energy − E_FCI (mHa)")
    axs[0].set_title("QSCI error does not follow the variational error", loc="left")
    axs[1].set_title("it follows the size of the sampled subspace", loc="left")
    proto = {"lin": "10^5 samples, 10 batches of 4,000, max_dim 4,000, one SQD iteration",
             "n2631g": "10^6 samples, one batch, no max_dim"}[a.proto]
    fig.text(0.01, 0.012, f"{nmol} molecules, norb 15-16, STO-3G. Lin et al. protocol: {proto}; QSCI = mean of the "
             f"independent draws.\nGray: one molecule and candidate. Blue: median per candidate. Log axes.",
             color=MUTED, fontsize=7.5, va="bottom")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, facecolor=SURF)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
