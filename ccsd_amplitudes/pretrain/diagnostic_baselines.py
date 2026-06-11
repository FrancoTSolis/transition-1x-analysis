#!/usr/bin/env python3
"""Diagnostic baselines: can ANY model predict kappa from t2?

Tests:
1. Trivial baseline: predict all zeros (MSE = Var(target))
2. Per-orbital-pair mean: predict mean kappa per (norb, pair_type)
3. Linear regression: flatten t2-slice → linear → kappa entry
4. Small MLP on t2-slices
5. Input/output feature scale analysis

This answers: is there signal in t2→kappa, or is the task fundamentally
too hard for the features we have?
"""

import os
import sys
import numpy as np
from pathlib import Path
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, mean_squared_error

def load_sample(job_dir):
    from pretrain.dataset import _parse_amplitude_file
    job_dir = Path(job_dir)
    t2 = _parse_amplitude_file(job_dir / "ccsd_t2.dat").astype(np.float32)
    kr = np.load(job_dir / "kappa_real.npy").astype(np.float32)
    ki = np.load(job_dir / "kappa_imag.npy").astype(np.float32)
    J = np.load(job_dir / "lucj_diag_coulomb_mats.npy").astype(np.float32)
    nocc, _, nvirt, _ = t2.shape
    return t2, kr[0], ki[0], J[0], nocc, nvirt


def main():
    jobs_dir = Path("jobs")
    job_dirs = sorted(
        d for d in jobs_dir.iterdir()
        if d.is_dir() and (d / "kappa_real.npy").exists()
    )

    # Use a fixed norb subset for fair comparison (most common size)
    # Find most common norb
    norb_counts = {}
    for d in job_dirs[:3000]:
        t2 = np.load(d / "kappa_real.npy")
        norb = t2.shape[1]
        norb_counts[norb] = norb_counts.get(norb, 0) + 1

    target_norb = max(norb_counts, key=norb_counts.get)
    print(f"Most common norb: {target_norb} ({norb_counts[target_norb]} samples)")
    print(f"Top 5 norbs: {sorted(norb_counts.items(), key=lambda x: -x[1])[:5]}")

    # Collect fixed-norb samples
    samples = []
    for d in job_dirs:
        kr = np.load(d / "kappa_real.npy")
        if kr.shape[1] != target_norb:
            continue
        try:
            t2, kr0, ki0, J0, nocc, nvirt = load_sample(d)
        except:
            continue
        samples.append((t2, kr0, ki0, J0, nocc, nvirt, d.name))
        if len(samples) >= 800:
            break

    n_samples = len(samples)
    print(f"\nCollected {n_samples} samples with norb={target_norb}")

    nocc = samples[0][4]
    nvirt = samples[0][5]
    norb = nocc + nvirt
    print(f"nocc={nocc}, nvirt={nvirt}, norb={norb}")

    # ============================================================
    # Feature scale analysis
    # ============================================================
    print("\n" + "=" * 60)
    print("INPUT/OUTPUT SCALE ANALYSIS")
    print("=" * 60)

    t2_vals = np.array([s[0].ravel() for s in samples[:200]])
    kr_vals = np.array([s[1].ravel() for s in samples[:200]])
    ki_vals = np.array([s[2].ravel() for s in samples[:200]])
    j_vals = np.array([s[3].ravel() for s in samples[:200]])

    print(f"t2:         mean={t2_vals.mean():.6f}, std={t2_vals.std():.6f}, range=[{t2_vals.min():.4f}, {t2_vals.max():.4f}]")
    print(f"kappa_real: mean={kr_vals.mean():.6f}, std={kr_vals.std():.6f}, range=[{kr_vals.min():.4f}, {kr_vals.max():.4f}]")
    print(f"kappa_imag: mean={ki_vals.mean():.6f}, std={ki_vals.std():.6f}, range=[{ki_vals.min():.4f}, {ki_vals.max():.4f}]")
    print(f"J:          mean={j_vals.mean():.6f}, std={j_vals.std():.6f}, range=[{j_vals.min():.4f}, {j_vals.max():.4f}]")

    # ============================================================
    # Baseline 0: Predict all zeros
    # ============================================================
    print("\n" + "=" * 60)
    print("BASELINE 0: Predict all zeros")
    print("=" * 60)

    all_kr = np.stack([s[1] for s in samples])  # (N, norb, norb)
    all_ki = np.stack([s[2] for s in samples])
    all_J = np.stack([s[3] for s in samples])

    mse_zero_kr = np.mean(all_kr ** 2)
    mse_zero_ki = np.mean(all_ki ** 2)
    mse_zero_J = np.mean(all_J ** 2)
    print(f"MSE(0, kappa_real) = {mse_zero_kr:.4f}")
    print(f"MSE(0, kappa_imag) = {mse_zero_ki:.4f}")
    print(f"MSE(0, J)          = {mse_zero_J:.4f}")
    print(f"MSE(0, kappa_total)= {mse_zero_kr + mse_zero_ki:.4f}")

    # ============================================================
    # Baseline 1: Per-entry mean across training set
    # ============================================================
    print("\n" + "=" * 60)
    print("BASELINE 1: Predict per-entry mean")
    print("=" * 60)

    n_train = int(0.8 * n_samples)
    train_kr = all_kr[:n_train]
    test_kr = all_kr[n_train:]
    train_ki = all_ki[:n_train]
    test_ki = all_ki[n_train:]

    mean_kr = train_kr.mean(axis=0)
    mean_ki = train_ki.mean(axis=0)
    mse_mean_kr = np.mean((test_kr - mean_kr) ** 2)
    mse_mean_ki = np.mean((test_ki - mean_ki) ** 2)
    print(f"MSE(mean, kappa_real) = {mse_mean_kr:.4f}  (vs zero: {mse_zero_kr:.4f})")
    print(f"MSE(mean, kappa_imag) = {mse_mean_ki:.4f}  (vs zero: {mse_zero_ki:.4f})")
    print(f"Per-entry mean helps: {mse_mean_kr < mse_zero_kr}")

    # ============================================================
    # Baseline 2: Linear regression t2_flat → kappa diagonal
    # ============================================================
    print("\n" + "=" * 60)
    print("BASELINE 2: Ridge regression on J diagonal (easiest target)")
    print("=" * 60)

    # J_aa off-diagonal (tridiagonal): just norb-1 values per sample
    all_j_offdiag = np.array([
        np.array([s[3][0, i, i+1] for i in range(norb-1)])
        for s in samples
    ])

    # Feature: flatten t2 (but it's huge). Use summary stats instead.
    all_t2_feats = np.array([
        np.concatenate([
            s[0].mean(axis=(2,3)).ravel(),  # mean over virt indices: (nocc, nocc)
            s[0].mean(axis=(0,1)).ravel(),  # mean over occ indices: (nvirt, nvirt)
            [s[0].std(), np.abs(s[0]).max(), np.linalg.norm(s[0].ravel())],
        ])
        for s in samples
    ])
    print(f"Feature dim: {all_t2_feats.shape[1]}")

    X_train = all_t2_feats[:n_train]
    X_test = all_t2_feats[n_train:]
    y_train = all_j_offdiag[:n_train]
    y_test = all_j_offdiag[n_train:]

    ridge = Ridge(alpha=1.0)
    ridge.fit(X_train, y_train)
    y_pred = ridge.predict(X_test)
    r2_j = r2_score(y_test, y_pred)
    mse_j = mean_squared_error(y_test, y_pred)
    print(f"J_aa off-diagonal: R² = {r2_j:.4f}, MSE = {mse_j:.6f}")
    print(f"  (zero baseline MSE = {np.mean(y_test**2):.6f})")

    # ============================================================
    # Baseline 3: Ridge regression on kappa diagonal
    # ============================================================
    print("\n" + "=" * 60)
    print("BASELINE 3: Ridge regression on kappa_real diagonal")
    print("=" * 60)

    all_kr_diag = np.array([np.diag(s[1]) for s in samples])  # should be ~0 (anti-sym)
    all_kr_offdiag1 = np.array([
        np.array([s[1][i, i+1] for i in range(norb-1)])
        for s in samples
    ])

    y_train_kr = all_kr_offdiag1[:n_train]
    y_test_kr = all_kr_offdiag1[n_train:]

    ridge_kr = Ridge(alpha=1.0)
    ridge_kr.fit(X_train, y_train_kr)
    y_pred_kr = ridge_kr.predict(X_test)
    r2_kr = r2_score(y_test_kr, y_pred_kr)
    mse_kr = mean_squared_error(y_test_kr, y_pred_kr)
    print(f"kappa_real 1st off-diagonal: R² = {r2_kr:.4f}, MSE = {mse_kr:.6f}")
    print(f"  (zero baseline MSE = {np.mean(y_test_kr**2):.6f})")

    # ============================================================
    # Baseline 4: Ridge on kappa_real full upper triangle
    # ============================================================
    print("\n" + "=" * 60)
    print("BASELINE 4: Ridge on kappa_real upper triangle (all entries)")
    print("=" * 60)

    triu_idx = np.triu_indices(norb, k=1)
    all_kr_triu = np.array([s[1][triu_idx] for s in samples])

    y_train_full = all_kr_triu[:n_train]
    y_test_full = all_kr_triu[n_train:]

    ridge_full = Ridge(alpha=10.0)
    ridge_full.fit(X_train, y_train_full)
    y_pred_full = ridge_full.predict(X_test)
    r2_full = r2_score(y_test_full.ravel(), y_pred_full.ravel())
    mse_full = mean_squared_error(y_test_full, y_pred_full)
    print(f"kappa_real upper tri: R² = {r2_full:.4f}, MSE = {mse_full:.6f}")
    print(f"  (zero baseline MSE = {np.mean(y_test_full**2):.6f})")

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Target':<30} {'Zero MSE':>10} {'Ridge MSE':>10} {'R²':>8} {'Signal?':>8}")
    print("-" * 70)
    print(f"{'J_aa off-diagonal':<30} {np.mean(y_test**2):>10.6f} {mse_j:>10.6f} {r2_j:>8.4f} {'YES' if r2_j > 0.05 else 'WEAK'}")
    print(f"{'kappa_real 1st off-diag':<30} {np.mean(y_test_kr**2):>10.6f} {mse_kr:>10.6f} {r2_kr:>8.4f} {'YES' if r2_kr > 0.05 else 'WEAK'}")
    print(f"{'kappa_real upper tri':<30} {np.mean(y_test_full**2):>10.6f} {mse_full:>10.6f} {r2_full:>8.4f} {'YES' if r2_full > 0.05 else 'WEAK'}")


if __name__ == "__main__":
    main()
