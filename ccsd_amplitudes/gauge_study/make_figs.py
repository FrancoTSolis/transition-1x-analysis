#!/usr/bin/env python3
"""Figures for the eigenvector-learnability study report.

Reads gauge_study/figs/exp2_pairs.npz and exp3_cache_pairs.npz.
Run from ccsd_amplitudes/:  python3 -m gauge_study.make_figs
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGS = "gauge_study/figs"

TEAL, RED, DARK, PURPLE, GREY = "#2a9d8f", "#e76f51", "#264653", "#7b2d8b", "#8d99ae"


def fig_ratios():
    e2 = np.load(f"{FIGS}/exp2_pairs.npz")
    e3 = np.load(f"{FIGS}/exp3_cache_pairs.npz")

    def ratio2(key):
        d = e2[key]
        return d[e2["is_close"]].mean() / d[~e2["is_close"]].mean()

    def ratio3(key):
        d_in, d = e3["t2"], e3[key]
        k = max(5, len(d_in) // 10)
        idx = np.argsort(d_in)[:k]
        return d[idx].mean() / d.mean()

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))

    labels2 = [("kappa raw\n(deck metric)", ratio2("d_kraw"), RED),
               ("U raw", ratio2("d_Uraw"), RED),
               ("U orbit\n(phase+perm)", ratio2("d_Ualign"), TEAL),
               ("$\\lambda_0 v_0$\n(invariant)", ratio2("d_H"), DARK),
               ("$B$\n(invariant)", ratio2("d_B"), DARK),
               ("$w$ sorted\n(eigenvalues)", ratio2("d_w"), DARK),
               ("Z raw\n(J control)", ratio2("d_Zraw"), TEAL)]
    labels3 = [("kappa raw\n(deck: 0.98)", ratio3("kraw"), RED),
               ("U raw", ratio3("Uraw"), RED),
               ("U phase-\naligned", ratio3("Uphase"), TEAL),
               ("U true gauge\n(+rev+swap)", ratio3("Urestr"), TEAL),
               ("t2-hat\n(invariant)", ratio3("t2hat"), DARK),
               ("Z raw\n(deck: 0.60)", ratio3("Zraw"), TEAL)]

    for ax, labels, title, inp in [
        (axes[0], labels2,
         "Exact-DF labels (n29 group, 1410 molecules)", ratio2("d_t2")),
        (axes[1], labels3,
         "Compressed LUCJ labels (rhf_cache, 42 molecules)", ratio3("t2")),
    ]:
        names = [l[0] for l in labels]
        vals = [l[1] for l in labels]
        cols = [l[2] for l in labels]
        ax.bar(range(len(vals)), vals, color=cols)
        ax.axhline(1.0, color="k", lw=0.8, ls="--")
        ax.axhline(inp, color=GREY, lw=1.2, ls=":",
                   label=f"input $t_2$ contrast ({inp:.2f})")
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(names, fontsize=8)
        ax.set_ylabel("closest-decile / random-pair distance ratio")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(0, 1.1)
        ax.legend(fontsize=8, loc="lower right")
        for i, v in enumerate(vals):
            ax.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    fig.suptitle("Same molecules, same labels, different coordinates: the invariant content of U is as smooth as the 'learnable' J",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{FIGS}/fig_ratios.png", dpi=130)
    plt.close(fig)


def fig_scatter():
    e2 = np.load(f"{FIGS}/exp2_pairs.npz")
    sel = np.random.default_rng(0).choice(len(e2["d_t2"]), 12000, replace=False)
    x = e2["d_t2"][sel]
    panels = [
        ("kappa raw entries (deck's fig2)", e2["d_kraw"][sel], RED),
        ("U orbit distance (gauge-aligned)", e2["d_Ualign"][sel], TEAL),
        ("invariant $B=\\lambda_0 v_0 v_0^T$", e2["d_B"][sel], DARK),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True)
    for ax, (name, y, c) in zip(axes, panels):
        ax.scatter(x, y, s=3, alpha=0.25, color=c, rasterized=True)
        r = np.corrcoef(x, y)[0, 1]
        ax.set_title(f"{name}\ncorr = {r:+.2f}", fontsize=10)
        ax.set_xlabel(r"$\|\Delta t_2\|$ between molecules")
    axes[0].set_ylabel("label distance")
    fig.suptitle("Exact-DF labels, n29 group: pairwise distances vs input distance",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{FIGS}/fig_scatter.png", dpi=130)
    plt.close(fig)


def fig_learnability():
    # numbers from exp4 log (ridge on t2-PCA256, 1410 molecules, 80/20)
    entries = [
        ("kappa raw\n(deck's target)", 0.176, RED),
        ("kappa template-\naligned", 0.362, TEAL),
        ("B entries\n(invariant)", 0.367, DARK),
        ("$\\lambda_0 v_0$ sign-\ncanonical (invariant)", 0.609, DARK),
        ("Z upper-tri\n(control)", 0.956, TEAL),
        ("inner eigenvalues\n$w$ (invariant)", 0.979, DARK),
    ]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    names = [e[0] for e in entries]
    vals = [e[1] for e in entries]
    cols = [e[2] for e in entries]
    ax.bar(range(len(vals)), vals, color=cols)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("test $R^2$ (pooled), ridge on t2-PCA256")
    ax.set_title("The learnability ladder (identical features & learner, only the target's\n"
                 "coordinates change): eigenvalues $\\gg$ invariant eigenvector content $\\gg$ raw coordinates",
                 fontsize=10)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.015, f"{v:+.2f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/fig_learnability.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    fig_ratios()
    fig_scatter()
    fig_learnability()
    print("figures written to", FIGS)
