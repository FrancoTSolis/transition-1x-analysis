#!/usr/bin/env python3
"""Deck v3 results figure: retargeted run vs the kappa baseline.

Left: gauge-aware validation metrics of the invariant run per epoch.
Right: invariant-content error vs the old model's zero-information line,
       plus the old kappa-MSE flatline (normalized) for contrast.

Run from ccsd_amplitudes/:  python3 -m gauge_study.make_v3_fig
"""
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEAL, RED, DARK, PURPLE = "#2a9d8f", "#e76f51", "#264653", "#7b2d8b"


def parse(path, keys):
    rows = []
    pat = re.compile(r"Epoch (\d+)/\d+\s+(.*)")
    for line in open(path):
        m = pat.search(line)
        if not m:
            continue
        row = {"epoch": int(m.group(1))}
        for k in keys:
            km = re.search(rf"{k}=([-\d.e+]+)", m.group(2))
            if km:
                row[k] = float(km.group(1))
        rows.append(row)
    return rows


def main():
    new = parse("train_invariant.log",
                ["val_kappa", "val_Berr_median", "val_cos_median",
                 "val_lam_r2", "val_z"])
    old = parse("train_pretrain.log", ["val_kappa"])

    ep = [r["epoch"] for r in new]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))

    ax = axes[0]
    ax.plot(ep, [r["val_cos_median"] for r in new], "o-", ms=3, color=TEAL,
            label=r"median $|\langle \hat v_0, v_0\rangle|$")
    ax.plot(ep, [r["val_lam_r2"] for r in new], "s-", ms=3, color=PURPLE,
            label=r"$\lambda_0$ $R^2$")
    ax.plot(ep, [r["val_kappa"] for r in new], "^-", ms=3, color=DARK,
            label=r"sign-invariant $v_0$ loss")
    ax.set_xlabel("epoch"); ax.set_ylim(-0.05, 1.05)
    ax.set_title("Invariant retarget: the 'unlearnable' rotation content trains\n"
                 "(3,020 val molecules, all sizes)", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(ep, [r["val_Berr_median"] for r in new], "o-", ms=3, color=TEAL,
            label=r"v3 invariant: median $\|\hat B - B\|/\|B\|$")
    ax.axhline(1.0, color=RED, lw=2, ls="--",
               label="v2 baseline (Z,κ)→t̂2: 1.000 (zero invariant content)")
    oe = [r["epoch"] for r in old]
    ok = np.array([r["val_kappa"] for r in old])
    ax.plot(oe, ok / ok[0], ":", color=RED, alpha=0.7,
            label="v2 κ-MSE, normalized to epoch 1 (flat = never learned)")
    ax.set_xlabel("epoch"); ax.set_ylim(0, 1.5)
    ax.set_title("Invariant-content error: below the baseline's\n"
                 "zero-information line from epoch 3", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.tight_layout()
    for out in ("gauge_study/figs/fig_v3_results.png",
                "slides/fig_v3_results.png", "docs/fig_v3_results.png"):
        fig.savefig(out, dpi=130)
    print("saved fig_v3_results.png")


if __name__ == "__main__":
    main()
