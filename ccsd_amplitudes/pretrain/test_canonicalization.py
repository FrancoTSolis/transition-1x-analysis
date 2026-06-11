#!/usr/bin/env python3
"""Test whether gauge-canonicalizing U makes kappa learnable.

Two experiments:
1. DETERMINISM: run ffsim's compressed DF on the same t2 twice; compare U.
   If U differs, the optimizer itself is non-deterministic (gauge orbit).
2. CANONICALIZATION: sign-fix + sort columns of U (MAP-style), recompute
   kappa = logm(U_canon), and re-test linear learnability vs the raw kappa.

Usage:
    python -m pretrain.test_canonicalization --cache-dir rhf_cache
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import logm
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score


def canonicalize_U(U):
    """MAP-style gauge fix of a unitary U (columns are eigenvector-like).

    1. Phase-fix each column: largest-magnitude entry made positive real.
    2. Sort columns by the row-index of their largest-magnitude entry.
    """
    n = U.shape[0]
    Uc = U.copy()
    # Phase-fix columns
    for j in range(n):
        col = Uc[:, j]
        k = np.argmax(np.abs(col))
        phase = col[k] / np.abs(col[k])
        Uc[:, j] = col / phase
    # Sort columns by argmax row index (canonical order), tie-break by that value
    keys = []
    for j in range(n):
        k = int(np.argmax(np.abs(Uc[:, j])))
        keys.append((k, -np.abs(Uc[k, j])))
    order = sorted(range(n), key=lambda j: keys[j])
    return Uc[:, order]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", default="rhf_cache")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    files = sorted(cache_dir.glob("*.npz"))
    samples = [np.load(f, allow_pickle=True) for f in files]
    n = len(samples)
    nocc = int(samples[0]["nocc"]); nvirt = int(samples[0]["nvirt"]); norb = int(samples[0]["norb"])
    print(f"{n} samples, norb={norb}, nocc={nocc}, nvirt={nvirt}")

    # ===== Experiment 1: determinism of the optimizer =====
    print("\n=== 1. DETERMINISM: same t2 -> same U? ===")
    import ffsim
    s0 = samples[0]
    t1, t2 = s0["t1"], s0["t2"]
    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    pairs_ab = [(p, p) for p in range(norb)]
    U_runs = []
    for trial in range(2):
        op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
            t2=t2, t1=t1, n_reps=2,
            interaction_pairs=(pairs_aa, pairs_ab),
            optimize=True, options=dict(maxiter=30),
        )
        U_runs.append(op.orbital_rotations[0])
    diff = np.abs(U_runs[0] - U_runs[1]).max()
    print(f"  max|U_run1 - U_run2| = {diff:.2e}")
    print(f"  Optimizer deterministic: {diff < 1e-6}")

    # ===== Experiment 2: canonicalization learnability =====
    print("\n=== 2. CANONICALIZATION: does gauge-fixing U help kappa learnability? ===")
    n_train = int(0.8 * n)

    # Build raw and canonicalized kappa from the stored U (reconstruct U = expm(kappa))
    raw_kr, can_kr = [], []
    for s in samples:
        # reconstruct U from stored kappa
        kappa = s["kappa_real"][0] + 1j * s["kappa_imag"][0]
        from scipy.linalg import expm
        U = expm(kappa)
        raw_kr.append(kappa.real)
        Uc = canonicalize_U(U)
        kappa_c = logm(Uc)
        can_kr.append(kappa_c.real)
    raw_kr = np.array(raw_kr)
    can_kr = np.array(can_kr)

    print(f"  raw kappa_re   std across molecules (per entry, mean): {raw_kr.std(axis=0).mean():.4f}")
    print(f"  canon kappa_re std across molecules (per entry, mean): {can_kr.std(axis=0).mean():.4f}")

    # Features: t2 summary
    feats = np.array([
        np.concatenate([s["t2"].mean(axis=(2,3)).ravel(),
                        s["t2"].mean(axis=(0,1)).ravel(),
                        [s["t2"].std(), np.abs(s["t2"]).max()]]) for s in samples
    ])
    Xtr, Xte = feats[:n_train], feats[n_train:]
    triu = np.triu_indices(norb, k=1)

    for lbl, K in [("RAW kappa_re", raw_kr), ("CANON kappa_re", can_kr)]:
        y = K[:, triu[0], triu[1]]
        r = Ridge(10.0).fit(Xtr, y[:n_train])
        r2 = r2_score(y[n_train:].ravel(), r.predict(Xte).ravel())
        print(f"  {lbl} upper-tri: linear R²={r2:.4f}")

    # Local per-pair test on canonicalized
    print("\n  Local per-pair (t2[i,j,:,:] -> kappa[i,j]):")
    for lbl, K in [("RAW", raw_kr), ("CANON", can_kr)]:
        r2s = []
        for i in range(nocc):
            for j in range(i + 1, nocc):
                X = np.array([s["t2"][i, j].reshape(-1) for s in samples])
                y = K[:, i, j]
                r = Ridge(0.1).fit(X[:n_train], y[:n_train])
                r2s.append(r2_score(y[n_train:], r.predict(X[n_train:])))
        print(f"    {lbl}: mean R²={np.mean(r2s):.4f}, max={np.max(r2s):.4f}")


if __name__ == "__main__":
    main()
