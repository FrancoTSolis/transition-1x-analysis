#!/usr/bin/env python3
"""Experiment 6: end-to-end proof of concept for the retargeted pipeline.

Instead of regressing U/kappa (gauge-scrambled), predict the invariant
sufficient statistic of the n_reps=2 exact-DF label:

    lam0 (signed scalar)  +  v0 (unit vector, sign-canonicalized)

then reconstruct (U, Z) at inference with the SAME deterministic linear
algebra ffsim uses (one-body matrix -> quadratures -> eigh). Score with
metrics that respect the gauge:

    t2hat error   ||t2hat(pred) - t2hat(label)|| / ||t2hat(label)||
    B error       same for the rank-1 invariant
    U orbit dist  min over phases+perm of ||U_pred - U_label||, weighted
                  comparison H = U diag(w) U^dag as well

Baselines: predict the dataset-mean v0; random train molecule's label.
Learner: ridge on t2-PCA256 (same features as exp4; deliberately weak --
the point is the target/metric, not the model).

Run from ccsd_amplitudes/:  python3 -m gauge_study.exp6_retarget_poc
"""
import time

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from .common import (align_unitaries, exact_df_t2, load_n29, quadrature,
                     vec_to_onebody)

SEED = 0


def build_from_v(lam, v, nocc, nvirt):
    """Deterministic (Z, U, H) reconstruction from predicted (lam, v)."""
    v = v / np.linalg.norm(v)
    M = vec_to_onebody(v, nocc, nvirt)
    Zs, Us, Hs = [], [], []
    for sign, coeff in ((1, lam), (-1, -lam)):
        H = quadrature(M, sign)
        w, U = np.linalg.eigh(H)
        Zs.append(coeff * np.outer(w, w))
        Us.append(U)
        Hs.append(H)
    return np.array(Zs), np.array(Us), np.array(Hs)


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    rows = load_n29()
    n = len(rows)
    nocc, nvirt = rows[0]["nocc"], rows[0]["nvirt"]
    order = rng.permutation(n)
    rows = [rows[i] for i in order]
    ntr = int(0.8 * n)

    lam = np.empty(n)
    V = np.empty((n, nocc * nvirt))
    for i, r in enumerate(rows):
        d = exact_df_t2(r["t2"], n_reps=2)
        v0 = d["v0"]
        v0 = v0 * np.sign(v0[np.argmax(np.abs(v0))])
        lam[i] = d["lam0"]
        V[i] = v0

    X = np.stack([r["t2"].ravel() for r in rows]).astype(np.float32)
    pca = PCA(n_components=256, svd_solver="randomized", random_state=SEED)
    Xtr = pca.fit_transform(X[:ntr])
    Xte = pca.transform(X[ntr:])

    r_lam = Ridge(1.0).fit(Xtr, lam[:ntr])
    lam_hat = r_lam.predict(Xte)
    r_v = Ridge(1.0).fit(Xtr, V[:ntr])
    V_hat = r_v.predict(Xte)
    print(f"lam0 ridge R^2 = {r2_score(lam[ntr:], lam_hat):+.3f}")
    print(f"v0   ridge R^2 (pooled) = "
          f"{r2_score(V[ntr:].ravel(), V_hat.ravel()):+.3f}")
    # sign-invariant v0 accuracy: |<v_hat, v_true>| after normalization
    Vn = V_hat / np.linalg.norm(V_hat, axis=1, keepdims=True)
    ov = np.abs(np.einsum("ij,ij->i", Vn, V[ntr:]))
    print(f"|<v0_hat, v0>| : median={np.median(ov):.3f}  "
          f"q25={np.quantile(ov, .25):.3f}  q75={np.quantile(ov, .75):.3f}")

    # mean-v0 baseline overlap
    vbar = V[:ntr].mean(0)
    vbar /= np.linalg.norm(vbar)
    ov_bar = np.abs(V[ntr:] @ vbar)
    print(f"baseline (train-mean v0) overlap: median={np.median(ov_bar):.3f}")

    # ---------------- reconstruct and score ----------------
    dU_pred, dU_rand = [], []
    tb_pred, tb_mean, tb_rand = [], [], []
    for t, itest in enumerate(range(ntr, n)):
        r = rows[itest]
        Zs, Us, _ = build_from_v(lam_hat[t], V_hat[t], nocc, nvirt)
        # orbit distance on rep 0, weighted implicitly by comparing H
        _, _, d = align_unitaries(r["U"][0], Us[0])
        dU_pred.append(d / np.linalg.norm(r["U"][0]))
        jrand = int(rng.integers(0, ntr))
        _, _, dr = align_unitaries(r["U"][0], rows[jrand]["U"][0])
        dU_rand.append(dr / np.linalg.norm(r["U"][0]))

        # invariant content error (B-space, equals t2hat error for 2 reps)
        B_true = lam[itest] * np.outer(V[itest], V[itest])
        B_pred = lam_hat[t] * np.outer(Vn[t], Vn[t])
        B_mean = lam[:ntr].mean() * np.outer(vbar, vbar)
        B_rand = lam[jrand] * np.outer(V[jrand], V[jrand])
        nrm = np.linalg.norm(B_true)
        tb_pred.append(np.linalg.norm(B_pred - B_true) / nrm)
        tb_mean.append(np.linalg.norm(B_mean - B_true) / nrm)
        tb_rand.append(np.linalg.norm(B_rand - B_true) / nrm)

    print("\n--- invariant-content error ||B_hat - B|| / ||B|| (median) ---")
    print(f"  ridge->reconstruct : {np.median(tb_pred):.3f}")
    print(f"  train-mean baseline: {np.median(tb_mean):.3f}")
    print(f"  random-label       : {np.median(tb_rand):.3f}")

    print("\n--- U orbit distance (rel), rep 0 (median) ---")
    print(f"  ridge->eigh pipeline: {np.median(dU_pred):.3f}")
    print(f"  random-label        : {np.median(dU_rand):.3f}")
    print("  (recall: relabeling the SAME molecule with a different LAPACK "
        "gives ~0.24 aligned -- that is the floor set by degenerate columns)")
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
