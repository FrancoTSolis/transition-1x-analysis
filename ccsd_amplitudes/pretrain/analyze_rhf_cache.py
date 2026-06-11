#!/usr/bin/env python3
"""Rich learnability analysis on the RHF-CCSD cache.

Tests whether t1/t2 can predict kappa and J with correct RHF data, using:
1. Per-pair t2 slice -> single kappa/J entry (strictly local linear model)
2. Full t2 -> kappa/J via MLP (sklearn MLPRegressor)
3. Orbital-energy-denominator features
4. J as a positive control

Usage:
    python -m pretrain.analyze_rhf_cache --cache-dir rhf_cache
"""

import argparse
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", default="rhf_cache")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    files = sorted(cache_dir.glob("*.npz"))
    samples = [np.load(f) for f in files]
    n = len(samples)
    if n < 20:
        print(f"Only {n} cached samples — need >=20. Wait for cache to build.")
        return

    nocc = int(samples[0]["nocc"])
    nvirt = int(samples[0]["nvirt"])
    norb = int(samples[0]["norb"])
    print(f"{n} samples, norb={norb}, nocc={nocc}, nvirt={nvirt}")

    n_train = int(0.8 * n)

    # Target scale
    all_kr = np.stack([s["kappa_real"][0] for s in samples])
    all_J = np.stack([s["J"][0, 0] for s in samples])
    print(f"\nkappa_real: std={all_kr.std():.4f}, per-elem var={np.mean(all_kr**2):.4f}")
    print(f"J_aa:       std={all_J.std():.4f}, per-elem var={np.mean(all_J**2):.4f}")

    # ===== Test A: strictly local — t2[i,j,:,:] -> kappa_re[i,j] (occ-occ) =====
    print("\n=== A. Local linear: t2[i,j,:,:] -> kappa_re[i,j] (occ pairs) ===")
    r2s_k, r2s_j = [], []
    for i in range(nocc):
        for j in range(i + 1, nocc):
            X = np.array([s["t2"][i, j].reshape(-1) for s in samples])
            yk = np.array([s["kappa_real"][0, i, j] for s in samples])
            yj = np.array([s["J"][0, 0, i, j] for s in samples])
            rk = Ridge(0.1).fit(X[:n_train], yk[:n_train])
            r2s_k.append(r2_score(yk[n_train:], rk.predict(X[n_train:])))
            rj = Ridge(0.1).fit(X[:n_train], yj[:n_train])
            r2s_j.append(r2_score(yj[n_train:], rj.predict(X[n_train:])))
    print(f"  kappa_re[i,j]: mean R²={np.mean(r2s_k):.4f}, max={np.max(r2s_k):.4f}")
    print(f"  J_aa[i,j]:     mean R²={np.mean(r2s_j):.4f}, max={np.max(r2s_j):.4f}")

    # ===== Test B: local occ-virt — t2[i,:,a,:] -> kappa_re[i, nocc+a] =====
    print("\n=== B. Local linear: t2[i,:,a,:] -> kappa_re[i,nocc+a] (occ-virt) ===")
    r2s = []
    for i in range(nocc):
        for a in range(nvirt):
            X = np.array([s["t2"][i, :, a, :].reshape(-1) for s in samples])
            y = np.array([s["kappa_real"][0, i, nocc + a] for s in samples])
            r = Ridge(0.1).fit(X[:n_train], y[:n_train])
            r2s.append(r2_score(y[n_train:], r.predict(X[n_train:])))
    print(f"  kappa_re occ-virt: mean R²={np.mean(r2s):.4f}, max={np.max(r2s):.4f}")

    # ===== Test C: full t2 -> kappa via MLP =====
    print("\n=== C. MLP: full t2 -> kappa_re upper-tri / J_aa off-diag ===")
    X_full = np.array([s["t2"].reshape(-1) for s in samples])
    scaler = StandardScaler().fit(X_full[:n_train])
    Xs = scaler.transform(X_full)
    Xtr, Xte = Xs[:n_train], Xs[n_train:]

    triu = np.triu_indices(norb, k=1)
    y_k = np.array([s["kappa_real"][0][triu] for s in samples])
    y_j = np.array([[s["J"][0, 0, i, i + 1] for i in range(norb - 1)] for s in samples])

    for lbl, y in [("kappa_re upper-tri", y_k), ("J_aa off-diag", y_j)]:
        mlp = MLPRegressor(hidden_layer_sizes=(256, 128), max_iter=500,
                           early_stopping=True, random_state=0)
        mlp.fit(Xtr, y[:n_train])
        pred = mlp.predict(Xte)
        r2 = r2_score(y[n_train:].ravel(), pred.ravel())
        print(f"  {lbl}: MLP R²={r2:.4f}")

    # ===== Test D: orbital energy denominators -> kappa =====
    print("\n=== D. Energy denominators (eps_i - eps_a) -> kappa occ-virt ===")
    eps = samples[0]["orb_energies"]
    print(f"  orbital energies available: {eps.shape}")
    X_eps = []
    for s in samples:
        e = s["orb_energies"]
        eo = e[:nocc]
        ev = e[nocc:]
        denom = (eo[:, None] - ev[None, :]).reshape(-1)  # (nocc*nvirt,)
        X_eps.append(denom)
    X_eps = np.array(X_eps)
    y_ov = np.array([s["kappa_real"][0][:nocc, nocc:].reshape(-1) for s in samples])
    r = Ridge(1.0).fit(X_eps[:n_train], y_ov[:n_train])
    r2 = r2_score(y_ov[n_train:].ravel(), r.predict(X_eps[n_train:]).ravel())
    print(f"  energy denom -> kappa occ-virt: R²={r2:.4f}")

    print("\n=== INTERPRETATION ===")
    print("If J shows positive R² but kappa does not, the kappa target itself")
    print("is the bottleneck (compressed-DF rotations are gauge-dependent / non-local).")


if __name__ == "__main__":
    main()
