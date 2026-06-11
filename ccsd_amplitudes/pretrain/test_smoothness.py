#!/usr/bin/env python3
"""Is kappa a smooth (Lipschitz) function of t2, like J is?

For all molecule pairs, correlate ||Δt2|| with ||Δkappa|| and ||ΔJ||.
A smooth target -> strong positive correlation (nearby inputs, nearby outputs).
A discontinuous target -> weak/no correlation.

Usage:
    python -m pretrain.test_smoothness --cache-dir rhf_cache
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", default="rhf_cache")
    args = parser.parse_args()

    files = sorted(Path(args.cache_dir).glob("*.npz"))
    S = [np.load(f) for f in files]
    n = len(S)
    nocc = int(S[0]["nocc"])

    t2 = [s["t2"].ravel() for s in S]
    kr = [s["kappa_real"][0].ravel() for s in S]
    Jaa = [s["J"][0, 0].ravel() for s in S]
    # reconstructed t2 from (U, Z) would be ideal; use J off-diag as gauge-invariant proxy

    d_t2, d_k, d_j = [], [], []
    for a, b in combinations(range(n), 2):
        d_t2.append(np.linalg.norm(t2[a] - t2[b]))
        d_k.append(np.linalg.norm(kr[a] - kr[b]))
        d_j.append(np.linalg.norm(Jaa[a] - Jaa[b]))
    d_t2 = np.array(d_t2); d_k = np.array(d_k); d_j = np.array(d_j)

    ck = np.corrcoef(d_t2, d_k)[0, 1]
    cj = np.corrcoef(d_t2, d_j)[0, 1]
    print(f"{n} molecules, {len(d_t2)} pairs\n")
    print(f"corr(||Δt2||, ||Δkappa||) = {ck:+.3f}   (smooth target -> high)")
    print(f"corr(||Δt2||, ||ΔJ_aa||)  = {cj:+.3f}   (smooth target -> high)")

    # Nearest-neighbor Lipschitz ratio: for the 10% closest-in-t2 pairs,
    # how big is the kappa jump vs J jump (normalized by their global std)?
    k_close = max(5, len(d_t2) // 10)
    idx = np.argsort(d_t2)[:k_close]
    print(f"\nAmong {k_close} closest-in-t2 pairs:")
    print(f"  mean ||Δt2||    = {d_t2[idx].mean():.4f}  (vs global {d_t2.mean():.4f})")
    print(f"  mean ||Δkappa|| = {d_k[idx].mean():.4f}  (vs global {d_k.mean():.4f})  "
          f"ratio={d_k[idx].mean()/d_k.mean():.2f}")
    print(f"  mean ||ΔJ_aa||  = {d_j[idx].mean():.4f}  (vs global {d_j.mean():.4f})  "
          f"ratio={d_j[idx].mean()/d_j.mean():.2f}")
    print("\nIf kappa ratio ~1 (no shrink) while J ratio <<1, kappa is NOT smooth in t2:")
    print("nearby t2 still gives far-apart kappa => discontinuous/gauge-entangled target.")


if __name__ == "__main__":
    main()
