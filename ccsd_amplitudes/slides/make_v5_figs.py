#!/usr/bin/env python3
"""Figures for deck_v5 (numbers from docs/chem_frame_pretraining.md, docs/rl_larger_molecules.md, docs/followups/)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parent
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 13, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF})


def style(ax):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)


# --- 1. pretraining: t2 compression residual on held-out n29 molecules (val141 median) -------------------------
rows = [("old MO-frame net, 1 call", 0.787), ("slot model, canonical frame, 1 call", 0.562),
        ("slot model, chemistry frame, 1 call", 0.367), ("slot model, chemistry frame, 4 calls", 0.331),
        ("slot model, chemistry frame, 8 calls", 0.317)]
label = 0.334
fig, ax = plt.subplots(figsize=(10, 4.2))
y = range(len(rows))[::-1]
ax.barh(list(y), [v for _, v in rows], height=0.55, color=BLUE)
for yi, (name, v) in zip(y, rows):
    ax.text(max(v, label) + 0.012, yi, f"{v:.3f}", va="center", color=INK, fontsize=12)
ax.set_yticks(list(y), [n for n, _ in rows])
ax.axvline(label, color=ORANGE, lw=2, ls="--")
ax.text(label + 0.008, len(rows) - 0.45, f"optimize=True labels {label:.3f}", color=INK2, fontsize=12)
ax.set_xlim(0, 0.9)
ax.set_xlabel("t2 compression residual, held-out n29 molecules (lower is better)")
ax.xaxis.grid(True, color=GRID)
ax.set_axisbelow(True)
style(ax)
fig.tight_layout()
fig.savefig(OUT / "fig_v5_pretrain.png", dpi=160)

# --- 2. energies by size: label vs pretrained vs NN+RL (rl4n29f) -----------------------------------------------
groups = ["norb 15–16\n(9, exact)", "norb 17–18\n(22, exact)", "norb 19\n(12, exact, unseen)",
          "norb 29 test\n(24, MPS χ256, unseen)"]
series = [("optimize=True labels", ORANGE, [62.0, 61.7, 66.7, 60.4]),
          ("pretrained, 4 calls", AQUA, [59.4, 51.1, 58.1, 54.5]),
          ("pretrained + RL, 4 calls", BLUE, [70.3, 68.6, 68.9, 66.5])]
fig, ax = plt.subplots(figsize=(10, 4.4))
off = {0: -0.18, 1: 0.0, 2: 0.18}
for k, (name, col, vals) in enumerate(series):
    xs = [i + off[k] for i in range(len(groups))]
    ax.scatter(xs, vals, s=120, color=col, label=name, zorder=3, edgecolors=SURF, linewidths=2)
    for x, v in zip(xs, vals):
        ax.text(x + 0.05, v, f"{v:.1f}", va="center", color=INK, fontsize=10.5)
for i in range(len(groups)):
    ax.plot([i + off[0], i + off[2]], [series[0][2][i], series[2][2][i]], color=GRID, lw=2, zorder=1)
ax.set_xticks(range(len(groups)), groups)
ax.set_ylim(48, 73)
ax.set_xlim(-0.5, len(groups) - 0.4)
ax.set_ylabel("% of CCSD correlation energy")
ax.yaxis.grid(True, color=GRID)
ax.set_axisbelow(True)
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.13), fontsize=12)
style(ax)
fig.tight_layout()
fig.savefig(OUT / "fig_v5_sizes.png", dpi=160)
print("ok")
