#!/usr/bin/env python3
"""Experiment 2: is the U target smooth in t2 once gauge is handled?

The deck's fig2 measured RAW ||delta kappa|| between molecules and found the
closest-in-t2 pairs are no closer in kappa than random pairs (ratio 0.98) --
and concluded U is unlearnable. Here we recompute that ratio with distances
that respect the label's symmetry group, on the n29 uniform-size set
(nocc=17, nvirt=12, exact-DF targets):

  raw          ||U_a - U_b||                     (what the deck measured)
  phase        min over per-column phases
  phase+perm   min over phases and column permutation (orbit distance)
  H            ||H_a - H_b||, H = U diag(w) U^T   (fully gauge-invariant)
  B            ||B_a - B_b||, B = lam0 v0 v0^T    (full invariant content)
  Z            raw + permutation-aligned Z (the "J is learnable" control)

Also: outer/inner spectral-gap statistics that govern where the eigenvector
map can be discontinuous (Davis-Kahan), and ratios conditioned on the
pair's spectral stability.

Run from ccsd_amplitudes/:  python3 -m gauge_study.exp2_smoothness
"""
import argparse
import time
from pathlib import Path

import numpy as np

from .common import (align_unitaries, exact_df_t2, load_n29,
                     nearest_decile_ratio, permute_Z)

FIGS = Path("gauge_study/figs")


def pairwise_dists(X):
    """Full pairwise Euclidean distance matrix via the Gram trick."""
    sq = np.einsum("ij,ij->i", X, X)
    D2 = sq[:, None] + sq[None, :] - 2.0 * (X @ X.T)
    np.fill_diagonal(D2, 0.0)
    return np.sqrt(np.maximum(D2, 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-mols", type=int, default=2100)
    ap.add_argument("--n-baseline", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    FIGS.mkdir(exist_ok=True)

    t0 = time.time()
    rows = load_n29(max_mols=args.max_mols)
    n = len(rows)
    nocc, nvirt, norb = rows[0]["nocc"], rows[0]["nvirt"], rows[0]["norb"]
    print(f"{n} molecules loaded ({time.time()-t0:.0f}s), "
          f"nocc={nocc} nvirt={nvirt}")

    # ---------------- per-molecule spectra + invariant content -------------
    lam0 = np.empty(n)
    gap_outer = np.empty(n)          # (|lam0|-|lam1|)/|lam0|
    V0 = np.empty((n, nocc * nvirt))
    n_zero_modes = np.empty(n, dtype=int)
    min_sigma_gap = np.empty(n)      # min adjacent gap of nonzero |w|, relative
    W = np.empty((n, norb))
    for i, r in enumerate(rows):
        d = exact_df_t2(r["t2"], n_reps=2)
        lam0[i] = d["lam0"]
        ae = np.abs(d["outer_eigs"])
        gap_outer[i] = (ae[0] - ae[1]) / ae[0]
        V0[i] = d["v0"]
        w = d["inner_eigs"][0]
        W[i] = np.sort(w)
        aw = np.sort(np.abs(w))[::-1]
        # spectrum is +/- sigma pairs + (nocc-nvirt) exact zeros; the basis
        # ambiguity of U is governed by gaps between DISTINCT sigmas
        sig = np.sort(w[w > 1e-9 * aw[0]])[::-1]
        n_zero_modes[i] = int((np.abs(w) < 1e-9 * aw[0]).sum())
        gaps = np.abs(np.diff(sig)) / sig[0]
        min_sigma_gap[i] = gaps.min()

    print("\n--- spectral statistics (n29, exact-DF labels) ---")
    print(f"outer relative gap (|lam0|-|lam1|)/|lam0|: "
          f"median={np.median(gap_outer):.3f}  "
          f"q10={np.quantile(gap_outer, 0.1):.3f}  "
          f"frac<0.05: {(gap_outer < 0.05).mean():.3f}  "
          f"frac<0.01: {(gap_outer < 0.01).mean():.3f}")
    print(f"inner zero modes per molecule: expected {nocc-nvirt}, "
          f"observed median={np.median(n_zero_modes):.0f} "
          f"(min={n_zero_modes.min()}, max={n_zero_modes.max()})")
    print(f"inner min adjacent sigma-gap (rel): median={np.median(min_sigma_gap):.4f}  "
          f"frac<0.01: {(min_sigma_gap < 0.01).mean():.3f}")

    # ---------------- pair selection by d_t2 --------------------------------
    X = np.stack([r["t2"].ravel() for r in rows]).astype(np.float64)
    D = pairwise_dists(X)
    iu = np.triu_indices(n, k=1)
    d_t2_all = D[iu]
    n_pairs = len(d_t2_all)

    k_close = int(0.1 * n_pairs)
    close_idx = np.argsort(d_t2_all)[:k_close]
    # same-reaction pairs (R/TS/P of one reaction): chemically-nearest pairs
    rxn = {}
    for i, r in enumerate(rows):
        rxn.setdefault(r["name"].rsplit("_", 1)[0], []).append(i)
    def lin(a, b):
        a, b = min(a, b), max(a, b)
        return a * n - a * (a + 1) // 2 + (b - a - 1)
    rxn_idx = np.array(sorted(
        lin(a, b) for members in rxn.values() if len(members) > 1
        for x, a in enumerate(members) for b in members[x + 1:]
    ), dtype=int)
    # cap the expensive aligned computations
    close_sel = (close_idx if k_close <= 40000
                 else rng.choice(close_idx, 40000, replace=False))
    base_sel = rng.choice(n_pairs, min(args.n_baseline, n_pairs), replace=False)
    sel = np.unique(np.concatenate([close_sel, base_sel, rxn_idx]))
    pairs = np.stack([iu[0][sel], iu[1][sel]], axis=1)
    is_close = np.isin(sel, close_idx)
    is_rxn = np.isin(sel, rxn_idx)
    print(f"\n{n_pairs} total pairs; aligning {len(sel)} "
          f"({is_close.sum()} in closest decile, {len(base_sel)} baseline, "
          f"{is_rxn.sum()} same-reaction)")

    # closed-form invariant distances on ALL selected pairs
    lam_a, lam_b = lam0[pairs[:, 0]], lam0[pairs[:, 1]]
    vdot = np.einsum("ij,ij->i", V0[pairs[:, 0]], V0[pairs[:, 1]])
    d_B = np.sqrt(np.maximum(lam_a**2 + lam_b**2 - 2 * lam_a * lam_b * vdot**2, 0))
    d_H = np.sqrt(np.maximum(
        lam_a**2 + lam_b**2 - 2 * np.abs(lam_a * lam_b * vdot), 0))
    d_w = np.linalg.norm(W[pairs[:, 0]] - W[pairs[:, 1]], axis=1)

    # loop for gauge-aligned distances
    d_t2 = d_t2_all[sel]
    m = len(pairs)
    d_kraw = np.empty(m); d_Uraw = np.empty(m)
    d_Uphase = np.empty(m); d_Ualign = np.empty(m)
    d_Zraw = np.empty(m); d_Zalign = np.empty(m)
    t0 = time.time()
    for t, (a, b) in enumerate(pairs):
        ra, rb = rows[a], rows[b]
        d_kraw[t] = np.linalg.norm(ra["kappa"] - rb["kappa"])
        d_Uraw[t] = np.linalg.norm(ra["U"][0] - rb["U"][0]) * np.sqrt(2)
        _, _, dp = align_unitaries(ra["U"][0], rb["U"][0], permute=False)
        d_Uphase[t] = dp * np.sqrt(2)
        _, perm, da = align_unitaries(ra["U"][0], rb["U"][0], permute=True)
        d_Ualign[t] = da * np.sqrt(2)
        d_Zraw[t] = np.linalg.norm(ra["Z"][0] - rb["Z"][0]) * np.sqrt(2)
        d_Zalign[t] = np.linalg.norm(
            ra["Z"][0] - permute_Z(rb["Z"][0], perm)) * np.sqrt(2)
        if (t + 1) % 20000 == 0:
            print(f"  aligned {t+1}/{m} pairs ({time.time()-t0:.0f}s)", flush=True)

    # ---------------- ratios ------------------------------------------------
    base = ~is_close
    def ratio(dv):
        return dv[is_close].mean() / dv[base].mean()

    print("\n--- closest-decile / random-pair ratio (lower = smoother) ---")
    names = [
        ("t2 (input, sanity)", d_t2),
        ("kappa RAW (deck's metric)", d_kraw),
        ("U RAW", d_Uraw),
        ("U phase-aligned", d_Uphase),
        ("U phase+perm aligned (orbit)", d_Ualign),
        ("lam0*v0, sign-aligned (invariant)", d_H),
        ("B = lam0 v0 v0^T  (full invariant)", d_B),
        ("inner eigenvalues w (sorted)", d_w),
        ("Z RAW (deck's J control)", d_Zraw),
        ("Z perm-aligned", d_Zalign),
    ]
    for name, dv in names:
        c = np.corrcoef(d_t2, dv)[0, 1]
        print(f"  {name:38s} ratio={ratio(dv):.3f}   corr(d_t2, .)={c:+.3f}")

    # same-reaction pairs: the deck's "near-identical t2" case by construction
    if is_rxn.sum():
        print(f"\n--- same-reaction pairs (n={is_rxn.sum()}): "
              f"mean distance / random-pair mean ---")
        for name, dv in names:
            print(f"  {name:38s} {dv[is_rxn].mean() / dv[base].mean():.3f}")

    # conditional on pair spectral stability
    stab = np.minimum(gap_outer[pairs[:, 0]], gap_outer[pairs[:, 1]])
    print("\n--- orbit/invariant ratio conditioned on outer-gap stability ---")
    for lo, hi in [(0.0, 0.05), (0.05, 0.2), (0.2, 1.0)]:
        msk = (stab >= lo) & (stab < hi)
        mc, mb = msk & is_close, msk & base
        if mc.sum() < 50 or mb.sum() < 50:
            continue
        print(f"  gap in [{lo},{hi}):  n_close={mc.sum():6d}  "
              f"U-orbit ratio={d_Ualign[mc].mean()/d_Ualign[mb].mean():.3f}  "
              f"B ratio={d_B[mc].mean()/d_B[mb].mean():.3f}  "
              f"kappa-raw ratio={d_kraw[mc].mean()/d_kraw[mb].mean():.3f}")

    np.savez(
        FIGS / "exp2_pairs.npz",
        d_t2=d_t2, d_kraw=d_kraw, d_Uraw=d_Uraw, d_Uphase=d_Uphase,
        d_Ualign=d_Ualign, d_B=d_B, d_H=d_H, d_w=d_w,
        d_Zraw=d_Zraw, d_Zalign=d_Zalign, is_close=is_close, stab=stab,
        gap_outer=gap_outer, min_sigma_gap=min_sigma_gap, lam0=lam0,
    )
    print(f"\nsaved pair data -> {FIGS/'exp2_pairs.npz'}  "
          f"({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
