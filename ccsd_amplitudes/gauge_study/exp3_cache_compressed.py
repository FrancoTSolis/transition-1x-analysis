#!/usr/bin/env python3
"""Experiment 3: decompose the non-smoothness of the COMPRESSED labels.

The deck's headline numbers (kappa ratio 0.98 vs J ratio 0.60) come from the
42-molecule rhf_cache whose labels are ffsim's optimize=True compressed DF
with LUCJ sparsity. Those labels have two independent pathologies:

  1. gauge: per-column phases always; permutations restricted by the LUCJ
     tridiagonal constraint to index reversal; rep swap.
  2. optimizer multimodality: L-BFGS from an (already gauge-arbitrary)
     exact-DF init, converging to input-sensitive local minima of a
     nonconvex objective that captures only ~6% of t2.

This experiment reproduces the deck's raw ratios, then progressively removes
gauge (phase; phase+reversal+swap; full Hungarian as an upper bound) and
finally compares the gauge-invariant physical content (reconstructed t2-hat).
Whatever non-smoothness remains after full alignment is optimizer noise, not
gauge -- it cannot be fixed by canonicalization, only by not regressing these
labels (reconstruction/energy losses, or predict-then-refine).

Run from ccsd_amplitudes/:  python3 -m gauge_study.exp3_cache_compressed
"""
from itertools import combinations

import numpy as np

from .common import (align_unitaries, exact_df_t2, load_cache,
                     nearest_decile_ratio, permute_Z, reconstruct_t2)


def align_restricted(Ua, Za, Ub, Zb):
    """min over per-column phases x {identity, reversal} x rep swap
    (the actual gauge group of the LUCJ-constrained compressed labels)."""
    n = Ua.shape[1]
    best = None
    rev = np.arange(n)[::-1]
    for order in ([0, 1], [1, 0]):
        for perm in (np.arange(n), rev):
            dU2 = dZ2 = 0.0
            for r, rb in enumerate(order):
                _, _, d = align_unitaries(Ua[r], Ub[rb][:, perm], permute=False)
                dU2 += d ** 2
                dZ2 += np.linalg.norm(Za[r] - permute_Z(Zb[rb], perm)) ** 2
            cand = (np.sqrt(dU2), np.sqrt(dZ2))
            if best is None or cand[0] < best[0]:
                best = cand
    return best


def main():
    rows = load_cache("rhf_cache", norb=22)
    # keep the dominant (nocc, nvirt) split
    from collections import Counter
    key = Counter((r["nocc"], r["nvirt"]) for r in rows).most_common(1)[0][0]
    rows = [r for r in rows if (r["nocc"], r["nvirt"]) == key]
    n = len(rows)
    nocc, nvirt = key
    print(f"{n} molecules (norb=22, nocc={nocc}, nvirt={nvirt}), "
          f"compressed LUCJ labels\n")

    # precompute per molecule: reconstruction t2hat and exact-DF init distance
    t2hat, d_init = [], []
    for r in rows:
        t2hat.append(reconstruct_t2(r["Z"], r["U"], nocc).real)
        init = exact_df_t2(r["t2"], n_reps=2)
        dU2 = 0.0
        for k in range(2):
            _, _, d = align_unitaries(init["U"][k], r["U"][k])
            dU2 += d ** 2
        d_init.append(np.sqrt(dU2) / np.linalg.norm(init["U"]))
        r["resid"] = (np.linalg.norm(t2hat[-1] - r["t2"])
                      / np.linalg.norm(r["t2"]))
    print(f"optimizer moved U from exact-DF init by (aligned, rel): "
          f"median={np.median(d_init):.3f}")
    print(f"reconstruction residual of stored labels: "
          f"median={np.median([r['resid'] for r in rows]):.3f}\n")

    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    cols = {k: np.empty(m) for k in
            ["t2", "kraw", "Uraw", "Uphase", "Urestr", "Uhung",
             "Zraw", "Zalign", "t2hat"]}
    for t, (a, b) in enumerate(pairs):
        ra, rb = rows[a], rows[b]
        cols["t2"][t] = np.linalg.norm(ra["t2"] - rb["t2"])
        cols["kraw"][t] = np.linalg.norm(ra["kappa"] - rb["kappa"])
        cols["Uraw"][t] = np.linalg.norm(ra["U"] - rb["U"])
        dp2 = 0.0
        for k in range(2):
            _, _, d = align_unitaries(ra["U"][k], rb["U"][k], permute=False)
            dp2 += d ** 2
        cols["Uphase"][t] = np.sqrt(dp2)
        dU, dZ = align_restricted(ra["U"], ra["Z"], rb["U"], rb["Z"])
        cols["Urestr"][t] = dU
        cols["Zalign"][t] = dZ
        dh2 = 0.0
        for k in range(2):
            _, _, d = align_unitaries(ra["U"][k], rb["U"][k], permute=True)
            dh2 += d ** 2
        cols["Uhung"][t] = np.sqrt(dh2)
        cols["Zraw"][t] = np.linalg.norm(ra["Z"] - rb["Z"])
        cols["t2hat"][t] = np.linalg.norm(t2hat[a] - t2hat[b])

    print("--- closest-decile / all-pair ratio (lower = smoother) ---")
    labels = [
        ("kappa RAW  (deck: 0.98)", "kraw"),
        ("U RAW", "Uraw"),
        ("U phase-aligned", "Uphase"),
        ("U phase+reversal+swap (true gauge)", "Urestr"),
        ("U + Hungarian perm (upper bound)", "Uhung"),
        ("Z RAW  (deck's J: 0.60)", "Zraw"),
        ("Z aligned", "Zalign"),
        ("t2hat  (gauge-invariant content)", "t2hat"),
    ]
    for name, k in labels:
        r = nearest_decile_ratio(cols["t2"], cols[k])
        c = np.corrcoef(cols["t2"], cols[k])[0, 1]
        print(f"  {name:38s} ratio={r:.3f}   corr={c:+.3f}")

    np.savez("gauge_study/figs/exp3_cache_pairs.npz", **cols)


if __name__ == "__main__":
    main()
