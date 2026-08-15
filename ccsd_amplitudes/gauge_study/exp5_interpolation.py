#!/usr/bin/env python3
"""Experiment 5: continuity of the label map along a t2 path.

Linearly interpolate t2 between two similar molecules (R and TS of one
reaction) and recompute the exact-DF labels at every step, exactly as the
label generator would. Tracks, between consecutive steps:

  raw ||dU||        : what supervised regression sees (LAPACK gauge included)
  orbit ||dU||      : after phase+permutation alignment
  ||dB||            : invariant content lam0 v0 v0^T
  outer gap         : (|lam0|-|lam1|)/|lam0|  (Davis-Kahan stability)

Expected picture: raw dU is erratic everywhere (gauge is redrawn at every
input); orbit dU is small and smooth except where the outer gap pinches
(genuine eigenvector crossing); dB is smooth wherever the gap is open.

Run from ccsd_amplitudes/:  python3 -m gauge_study.exp5_interpolation
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .common import align_unitaries, exact_df_t2, load_n29

FIGS = "gauge_study/figs"


def path_scan(t2a, t2b, n_steps=240):
    ss = np.linspace(0.0, 1.0, n_steps)
    Us, Bs, gaps, lam0s = [], [], [], []
    for s in ss:
        t2 = (1 - s) * t2a + s * t2b
        d = exact_df_t2(t2, n_reps=2)
        ae = np.abs(d["outer_eigs"])
        gaps.append((ae[0] - ae[1]) / ae[0])
        lam0s.append(d["lam0"])
        Us.append(d["U"][0])
        Bs.append(d["lam0"] * np.outer(d["v0"], d["v0"]))
    dU_raw, dU_orb, dB = [], [], []
    for t in range(n_steps - 1):
        dU_raw.append(np.linalg.norm(Us[t + 1] - Us[t]))
        _, _, d = align_unitaries(Us[t], Us[t + 1])
        dU_orb.append(d)
        dB.append(np.linalg.norm(Bs[t + 1] - Bs[t]))
    return ss, np.array(gaps), np.array(dU_raw), np.array(dU_orb), np.array(dB)


def main():
    rows = load_n29(with_U=False)
    by_name = {r["name"]: r for r in rows}

    # pick a same-reaction R->TS pair, plus the globally least-stable molecule
    # paired with its nearest neighbor (crossing likely on that segment)
    names = sorted(by_name)
    pair1 = None
    for nm in names:
        if nm.endswith("_R") and nm[:-2] + "_TS" in by_name:
            pair1 = (nm, nm[:-2] + "_TS")
            break

    gaps = {nm: None for nm in names}
    lam_all = []
    for nm in names:
        d = exact_df_t2(by_name[nm]["t2"], n_reps=2)
        ae = np.abs(d["outer_eigs"])
        gaps[nm] = (ae[0] - ae[1]) / ae[0]
    worst = min(names, key=lambda nm: gaps[nm])
    X = np.stack([by_name[nm]["t2"].ravel() for nm in names])
    xw = by_name[worst]["t2"].ravel()
    d2 = np.einsum("ij,ij->i", X - xw[None, :], X - xw[None, :])
    d2[names.index(worst)] = np.inf
    pair2 = (worst, names[int(np.argmin(d2))])
    print(f"pair1 (generic):  {pair1[0]} -> {pair1[1]}")
    print(f"pair2 (near-degenerate, gap={gaps[worst]:.4f}): "
          f"{pair2[0]} -> {pair2[1]}")

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex="col")
    for c, (na, nb) in enumerate([pair1, pair2]):
        ss, gap, dU_raw, dU_orb, dB = path_scan(
            by_name[na]["t2"], by_name[nb]["t2"])
        sm = 0.5 * (ss[1:] + ss[:-1])
        ax = axes[0, c]
        ax.semilogy(sm, dU_raw, lw=1, color="#e76f51",
                    label=r"raw $\|\Delta U\|$ (what MSE sees)")
        ax.semilogy(sm, dU_orb, lw=1, color="#2a9d8f",
                    label=r"orbit $\|\Delta U\|$ (gauge-aligned)")
        ax.semilogy(sm, dB, lw=1, color="#264653",
                    label=r"$\|\Delta B\|$ (invariant)")
        ax.set_title(f"{'generic pair' if c == 0 else 'near-degenerate pair'}"
                     f"\n{na} $\\to$ {nb}", fontsize=9)
        ax.set_ylabel("step-to-step change" if c == 0 else "")
        ax.legend(fontsize=7)
        ax2 = axes[1, c]
        ax2.semilogy(ss, gap, lw=1, color="#7b2d8b")
        ax2.set_ylabel(r"outer gap $(|\lambda_0|-|\lambda_1|)/|\lambda_0|$"
                       if c == 0 else "")
        ax2.set_xlabel("interpolation s")
        print(f"  {na}->{nb}: raw dU median {np.median(dU_raw):.3f}, "
              f"orbit dU median {np.median(dU_orb):.4f}, "
              f"dB median {np.median(dB):.4f}, min gap {gap.min():.4f}")
    fig.suptitle("Label continuity along a t2 path: raw labels jump everywhere; "
                 "the orbit/invariant is smooth away from spectral crossings",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/exp5_interpolation.png", dpi=130)
    print(f"saved {FIGS}/exp5_interpolation.png")


if __name__ == "__main__":
    main()
