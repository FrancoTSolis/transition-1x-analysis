#!/usr/bin/env python3
"""Experiment 4: learnability retest with gauge handled, n29 exact labels.

The deck's learnability test regressed RAW kappa entries and got R^2 = -0.16,
concluding "U is unlearnable". Here we regress, with identical features and
identical learners:

  (a) kappa RAW            -- reproduces the deck's negative result
  (b) kappa template-aligned  (all U's gauge-aligned to one train molecule)
  (c) lam0*v0 sign-canonicalized  (the minimal sufficient target: everything
      else in (U,Z) is deterministic linear algebra on top of it)
  (d) B = lam0 v0 v0^T entries    (sign-free invariant, sampled)
  (e) w  (inner eigenvalues, sorted -- spectrum of the quadrature)
  (f) Z entries            -- the deck's learnable control
plus
  (g) 1-NN label transfer scored by ORBIT distance vs raw distance
      (same predictor, two metrics -> shows the metric, not the signal,
       drives the "unlearnable" verdict).

Run from ccsd_amplitudes/:  python3 -m gauge_study.exp4_learnability
"""
import time

import numpy as np
from scipy.linalg import schur
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from .common import align_unitaries, exact_df_t2, load_n29

SEED = 0


def unitary_logm(U):
    """Principal log of a unitary matrix via complex Schur (robust, fast)."""
    T, Q = schur(U.astype(complex), output="complex")
    return Q @ np.diag(1j * np.angle(np.diag(T))) @ Q.conj().T


def kappa_feats(k):
    iu = np.triu_indices(k.shape[0], 1)
    iud = np.triu_indices(k.shape[0], 0)
    return np.concatenate([k.real[iu], k.imag[iud]])


def fit_ridge(Xtr, ytr, Xte, yte):
    best = -np.inf
    n_val = max(1, int(0.15 * len(Xtr)))
    for a in (0.1, 1.0, 10.0, 100.0, 1000.0):
        r = Ridge(a).fit(Xtr[:-n_val], ytr[:-n_val])
        v = r2_score(ytr[-n_val:].ravel(), r.predict(Xtr[-n_val:]).ravel())
        if v > best:
            best, keep = v, a
    r = Ridge(keep).fit(Xtr, ytr)
    return r2_score(yte.ravel(), r.predict(Xte).ravel()), keep


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    rows = load_n29()
    n = len(rows)
    nocc, nvirt, norb = rows[0]["nocc"], rows[0]["nvirt"], rows[0]["norb"]
    perm_mols = rng.permutation(n)
    rows = [rows[i] for i in perm_mols]
    ntr = int(0.8 * n)
    print(f"{n} molecules, train {ntr} / test {n-ntr}  ({time.time()-t0:.0f}s)")

    # ---------- per-molecule quantities ----------
    lamv, Bs, Ws, Zs = [], [], [], []
    for r in rows:
        d = exact_df_t2(r["t2"], n_reps=2)
        v0 = d["v0"]
        v0 = v0 * np.sign(v0[np.argmax(np.abs(v0))])       # sign canonicalize
        lamv.append(d["lam0"] * v0)
        B = d["lam0"] * np.outer(v0, v0)
        Bs.append(B)
        Ws.append(np.sort(d["inner_eigs"][0]))
        Zs.append(r["Z"][0])
    lamv = np.array(lamv)                                   # (n, 204)
    Ws = np.array(Ws)
    iu29 = np.triu_indices(norb, 0)
    Zt = np.array([z[iu29] for z in Zs])
    B_idx = rng.choice(nocc * nvirt * nocc * nvirt, 500, replace=False)
    Bt = np.array([b.ravel()[B_idx] for b in Bs])

    # kappa raw and template-aligned
    k_raw = np.array([kappa_feats(r["kappa"][0]) for r in rows])
    template = rows[0]["U"][0]
    k_ali_mats = []
    for r in rows:
        Ual, _, _ = align_unitaries(template, r["U"][0])
        k_ali_mats.append(unitary_logm(Ual))
    k_ali = np.array([kappa_feats(k) for k in k_ali_mats])
    print(f"targets built ({time.time()-t0:.0f}s)")

    # ---------- features: t2 PCA-256 ----------
    X = np.stack([r["t2"].ravel() for r in rows]).astype(np.float32)
    pca = PCA(n_components=256, svd_solver="randomized", random_state=SEED)
    Xtr_p = pca.fit_transform(X[:ntr])
    Xte_p = pca.transform(X[ntr:])
    print(f"PCA done, explained var {pca.explained_variance_ratio_.sum():.3f} "
          f"({time.time()-t0:.0f}s)")

    print("\n--- pooled test R^2, ridge on t2-PCA256 features ---")
    for name, Y in [
        ("kappa RAW (deck's target)", k_raw),
        ("kappa template-aligned", k_ali),
        ("lam0*v0 sign-canonical (204)", lamv),
        ("B entries (500 sampled)", Bt),
        ("inner eigenvalues w (29)", Ws),
        ("Z upper-tri (control)", Zt),
    ]:
        r2, a = fit_ridge(Xtr_p, Y[:ntr], Xte_p, Y[ntr:])
        print(f"  {name:34s} R^2 = {r2:+.3f}   (alpha={a})")

    # ---------- deck-style strictly local linear test ----------
    print("\n--- local linear test: t2[i,j,:,:] -> target[i,j] (occ pairs) ---")
    r2s_raw, r2s_ali, r2s_B = [], [], []
    for i in range(nocc):
        for j in range(i + 1, nocc):
            Xl = np.array([r["t2"][i, j].ravel() for r in rows])
            for tgt, out in [
                (np.array([r["kappa"][0].real[i, j] for r in rows]), r2s_raw),
                (np.array([k.real[i, j] for k in k_ali_mats]), r2s_ali),
                (np.array([b.reshape(nocc, nvirt, nocc, nvirt)[i, :, j, :].ravel()
                           for b in Bs]), r2s_B),
            ]:
                rr = Ridge(0.1).fit(Xl[:ntr], tgt[:ntr])
                out.append(r2_score(
                    np.asarray(tgt[ntr:]).ravel(),
                    rr.predict(Xl[ntr:]).ravel()))
    print(f"  kappa_re RAW      mean R^2 = {np.mean(r2s_raw):+.3f}")
    print(f"  kappa_re ALIGNED  mean R^2 = {np.mean(r2s_ali):+.3f}")
    print(f"  B[i,:,j,:] local  mean R^2 = {np.mean(r2s_B):+.3f}")

    # ---------- (g) 1-NN label transfer, two metrics ----------
    print("\n--- 1-NN label transfer (predict U of test = U of nearest train) ---")
    d2 = (np.einsum("ij,ij->i", X[ntr:], X[ntr:])[:, None]
          + np.einsum("ij,ij->i", X[:ntr], X[:ntr])[None, :]
          - 2.0 * X[ntr:] @ X[:ntr].T)
    nn = np.argmin(d2, axis=1)
    raw_nn, orbit_nn, raw_rnd, orbit_rnd = [], [], [], []
    for t, itest in enumerate(range(ntr, n)):
        Ut = rows[itest]["U"][0]
        for src, raw_out, orb_out in [
            (rows[int(nn[t])]["U"][0], raw_nn, orbit_nn),
            (rows[int(rng.integers(0, ntr))]["U"][0], raw_rnd, orbit_rnd),
        ]:
            raw_out.append(np.linalg.norm(Ut - src))
            _, _, d = align_unitaries(Ut, src)
            orb_out.append(d)
    print(f"  RAW   distance:  1-NN {np.mean(raw_nn):.3f}  vs random "
          f"{np.mean(raw_rnd):.3f}   (ratio {np.mean(raw_nn)/np.mean(raw_rnd):.3f})")
    print(f"  ORBIT distance:  1-NN {np.mean(orbit_nn):.3f}  vs random "
          f"{np.mean(orbit_rnd):.3f}   (ratio {np.mean(orbit_nn)/np.mean(orbit_rnd):.3f})")
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
