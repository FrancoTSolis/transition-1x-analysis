#!/usr/bin/env python3
"""Experiment 7: are canonical-gauge optimize=True (U, J) labels learnable?

Question (Sep 2026): the invariant retarget works because it removed the
gauge.  If the gauge of the *compressed* labels is fixed at the source
(canonical init before L-BFGS, generate_compressed_targets.py), is the map
t2 -> (U_opt, Z_opt) smooth and learnable, so that U and J can be regressed
directly in the optimize=True setting?

For one config directory <out-dir>/<config>/ on a same-shape molecule group:

  A. label facts: residual init vs opt, optimizer displacement |dkappa|,
     eigenphase range of U_init^dag U_opt (log branch-cut risk), and -- if the
     generator ran with --gauge-check -- the optimizer's gauge sensitivity
     (aligned distance between runs from differently-phased inits).
  B. smoothness: closest-decile / all-pairs distance ratio (lower = smoother;
     compare to the t2 input contrast) for
        U_opt raw | U_opt phase-aligned | H(U_opt) phase-free encoding |
        dkappa | Z_opt | dZ | t2hat_opt (invariant content) |
        U_init canonical raw / phase-aligned | B = lam0 v0 v0^T
  C. learnability: ridge on t2-PCA256, pooled test R^2 for
        kappa_opt (canonical) | dkappa | H(U_opt) | Z_opt (allowed entries) |
        dZ | t2hat_opt | lam0*v0 (invariant baseline) | Z_init
     plus 1-NN label transfer under raw vs phase-orbit distance.

Run from ccsd_amplitudes/:
    python3 -m gauge_study.exp7_compressed_canonical --config square_reg0
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from .common import nearest_decile_ratio
from .compressed_canonical import reconstruct_t2, unitary_log, z_mask

SEED = 0


# ------------------------------------------------------------- utilities

def pair_sqdist(X: np.ndarray) -> np.ndarray:
    """All-pairs squared Euclidean distances of rows (real or complex)."""
    Xf = X.reshape(len(X), -1)
    sq = np.einsum("ij,ij->i", Xf.conj(), Xf).real
    G = (Xf.conj() @ Xf.T).real
    D = sq[:, None] + sq[None, :] - 2 * G
    return np.maximum(D, 0.0)


def pair_phase_aligned_sqdist(U: np.ndarray, chunk: int = 200) -> np.ndarray:
    """All-pairs ||U_a - U_b D||^2 minimized over per-column phases D, summed
    over reps.  U: (N, reps, n, n).  ||u-v e^{i th}||^2 = 2 - 2|<u,v>|."""
    N, R, n, _ = U.shape
    D = np.zeros((N, N))
    for r in range(R):
        Ur = U[:, r]                                   # (N, n, n)
        for s in range(0, N, chunk):
            a = Ur[s:s + chunk]                        # (c, n, n)
            ov = np.einsum("aip,bip->abp", a.conj(), Ur)   # (c, N, n)
            D[s:s + chunk] += (2.0 - 2.0 * np.abs(ov)).sum(-1)
    return np.maximum(D, 0.0)


def upper(x: np.ndarray, k: int = 0) -> np.ndarray:
    iu = np.triu_indices(x.shape[-1], k)
    return x[..., iu[0], iu[1]]


def kappa_feats(k: np.ndarray) -> np.ndarray:
    """anti-Hermitian (reps,n,n) -> real features: real upper (k=1), imag upper (k=0)."""
    return np.concatenate([upper(k.real, 1).ravel(), upper(k.imag, 0).ravel()])


def herm_feats(h: np.ndarray) -> np.ndarray:
    return np.concatenate([upper(h.real, 0).ravel(), upper(h.imag, 1).ravel()])


def h_encoding(U: np.ndarray) -> np.ndarray:
    """Phase-free encoding of U keeping column order: H = U diag(0..n-1) U^dag."""
    n = U.shape[-1]
    w = np.arange(n, dtype=float) / (n - 1)
    return np.einsum("...ip,p,...jp->...ij", U, w, U.conj())


def fit_ridge(Xtr, ytr, Xte, yte):
    best, keep = -np.inf, 1.0
    n_val = max(1, int(0.15 * len(Xtr)))
    for a in (0.1, 1.0, 10.0, 100.0, 1000.0):
        r = Ridge(a).fit(Xtr[:-n_val], ytr[:-n_val])
        v = r2_score(ytr[-n_val:].ravel(), r.predict(Xtr[-n_val:]).ravel())
        if v > best:
            best, keep = v, a
    r = Ridge(keep).fit(Xtr, ytr)
    pred = r.predict(Xte)
    # pooled R^2 over all entries is inflated by between-entry mean differences
    # (a per-entry constant predictor scores high); also report the R^2 of the
    # molecule-specific part: both truth and prediction centered by the
    # per-entry TRAIN mean.
    mu = ytr.mean(0, keepdims=True)
    r2_pooled = r2_score(yte.ravel(), pred.ravel())
    ss_res = ((yte - pred) ** 2).sum()
    ss_tot = ((yte - mu) ** 2).sum()
    r2_centered = 1.0 - ss_res / max(ss_tot, 1e-30)
    return (r2_pooled, r2_centered), keep, pred


# ------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="rhf_targets_compressed")
    ap.add_argument("--config", default="square_reg0")
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--names-file", default="gauge_study/names_n29_16_13.txt")
    ap.add_argument("--max-mols", type=int, default=0)
    ap.add_argument("--save", default=None, help="npz path for pair columns")
    args = ap.parse_args()
    t0 = time.time()

    conn = args.config.rsplit("_reg", 1)[0]
    names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    cdir = Path(args.out_dir) / args.config
    idir = Path(args.out_dir) / "init"
    names = [n for n in names if (cdir / f"{n}.npz").exists() and (idir / f"{n}.npz").exists()]
    if args.max_mols:
        names = names[:args.max_mols]
    rng = np.random.default_rng(SEED)
    names = [names[i] for i in rng.permutation(len(names))]
    N = len(names)
    print(f"config={args.config}: {N} molecules with labels  ({time.time()-t0:.0f}s)")

    # ---------------- load ----------------
    t2s, Uo, Zo, Ui, Zi, dk, lamv, gaps, resid_o, resid_i = [], [], [], [], [], [], [], [], [], []
    nit, succ, gdist, gdist_raw, gdZ, t_opt = [], [], [], [], [], []
    for n in names:
        d = np.load(Path(args.data_dir) / f"{n}.npz")
        t2s.append(d["t2"].astype(np.float64))
        c = np.load(cdir / f"{n}.npz")
        i = np.load(idir / f"{n}.npz")
        Uo.append(c["U_re"] + 1j * c["U_im"]); Zo.append(c["Z"])
        Ui.append(i["U_re"] + 1j * i["U_im"]); Zi.append(i["Z"])
        dk.append(c["dkappa_re"] + 1j * c["dkappa_im"])
        lamv.append(float(i["lam0"]) * i["v0"].ravel()); gaps.append(float(i["gap"]))
        resid_o.append(float(c["resid"])); nit.append(int(c["nit"])); succ.append(bool(c["success"]))
        t_opt.append(float(c["time"]))
        resid_i.append(float(i["resid_init_square" if conn == "square" else "resid_init_all"]))
        if "gauge_dist" in c:
            gdist.append(float(c["gauge_dist"])); gdist_raw.append(float(c["gauge_dist_raw"]))
            gdZ.append(float(c["gauge_dZ"]))
    t2s = np.array(t2s); Uo = np.array(Uo); Zo = np.array(Zo); Ui = np.array(Ui); Zi = np.array(Zi)
    dk = np.array(dk); lamv = np.array(lamv); gaps = np.array(gaps)
    nocc, _, nvirt, _ = t2s.shape[1:]
    norb = nocc + nvirt
    mask = z_mask(conn, norb)
    Zi_m = Zi * mask[None, None]
    dZ = Zo - Zi_m
    t2hat = np.array([reconstruct_t2(Zo[i], Uo[i], nocc) for i in range(N)])
    B = np.einsum("ni,nj->nij", lamv, lamv) / np.maximum(
        np.linalg.norm(lamv, axis=1)[:, None, None] ** 2 / np.abs(lamv).max(1)[:, None, None] * 0 + 1, 1)
    # B = lam0 v0 v0^T with unit v0: lamv = lam0 v0 -> lam0 = ||lamv|| * sign
    lam0 = np.linalg.norm(lamv, axis=1)
    B = np.einsum("ni,nj->nij", lamv, lamv) / lam0[:, None, None]
    print(f"loaded ({time.time()-t0:.0f}s)  nocc={nocc} nvirt={nvirt}")

    # ---------------- A. label facts ----------------
    print("\n--- A. label facts ---")
    print(f"  residual  init median {np.median(resid_i):.3f}  ->  opt median {np.median(resid_o):.3f}"
          f"  (q25 {np.quantile(resid_o,.25):.3f}, q75 {np.quantile(resid_o,.75):.3f})")
    print(f"  L-BFGS: nit median {np.median(nit):.0f}  success frac {np.mean(succ):.2f}"
          f"  time/mol median {np.median(t_opt):.1f}s")
    dk_norm = np.linalg.norm(dk.reshape(N, -1), axis=1) / np.sqrt(2 * norb)
    eig_ph = np.array([np.max(np.abs(np.angle(np.linalg.eigvals(Ui[i, k].conj().T @ Uo[i, k]))))
                       for i in range(N) for k in range(2)]).reshape(N, 2).max(1)
    print(f"  |dkappa|_F/sqrt(2n): median {np.median(dk_norm):.3f}  q90 {np.quantile(dk_norm,.9):.3f}")
    print(f"  max eigenphase of U_init^dag U_opt: median {np.median(eig_ph):.2f} rad,"
          f"  frac > 0.9*pi: {np.mean(eig_ph > 0.9*np.pi):.3f}  (branch-cut risk of log)")
    print(f"  ||Z_opt||: median {np.median(np.linalg.norm(Zo.reshape(N,-1),axis=1)):.3f}"
          f"   ||Z_init(masked)||: {np.median(np.linalg.norm(Zi_m.reshape(N,-1),axis=1)):.3f}"
          f"   ||dZ||: {np.median(np.linalg.norm(dZ.reshape(N,-1),axis=1)):.3f}")
    if gdist:
        gd = np.array(gdist); gdr = np.array(gdist_raw); gz = np.array(gdZ)
        print(f"  optimizer gauge sensitivity (rerun from random-phase init):")
        print(f"    raw ||U-U'||: median {np.median(gdr):.3f}   phase-aligned: median {np.median(gd):.3f}"
              f"  q90 {np.quantile(gd,.9):.3f}   [scale ||U||=sqrt(2n)={np.sqrt(2*norb):.2f}]"
              f"   ||Z-Z'||: median {np.median(gz):.4f}")

    # ---------------- B. smoothness ----------------
    print("\n--- B. closest-decile / all-pairs ratio (lower = smoother) ---")
    iu = np.triu_indices(N, 1)
    d_in = np.sqrt(pair_sqdist(t2s))[iu]
    cols = {}
    cols["t2 (input)"] = d_in
    cols["U_init canonical RAW"] = np.sqrt(pair_sqdist(Ui))[iu]
    cols["U_init phase-aligned"] = np.sqrt(pair_phase_aligned_sqdist(Ui))[iu]
    cols["U_opt RAW"] = np.sqrt(pair_sqdist(Uo))[iu]
    cols["U_opt phase-aligned"] = np.sqrt(pair_phase_aligned_sqdist(Uo))[iu]
    cols["H(U_opt) phase-free"] = np.sqrt(pair_sqdist(h_encoding(Uo)))[iu]
    cols["dkappa (rel. to init)"] = np.sqrt(pair_sqdist(dk))[iu]
    cols["Z_opt"] = np.sqrt(pair_sqdist(Zo))[iu]
    cols["dZ = Z_opt - Z_init"] = np.sqrt(pair_sqdist(dZ))[iu]
    cols["Z_init (masked)"] = np.sqrt(pair_sqdist(Zi_m))[iu]
    cols["t2hat_opt (invariant content)"] = np.sqrt(pair_sqdist(t2hat))[iu]
    cols["B = lam0 v0 v0^T"] = np.sqrt(pair_sqdist(B))[iu]
    for k, v in cols.items():
        r = nearest_decile_ratio(d_in, v)
        c = np.corrcoef(d_in, v)[0, 1]
        print(f"  {k:32s} ratio={r:.3f}   corr={c:+.3f}")
    print(f"  ({time.time()-t0:.0f}s)")

    # ---------------- C. learnability ----------------
    print("\n--- C. pooled test R^2, ridge on t2-PCA256 (80/20 split) ---")
    ntr = int(0.8 * N)
    X = t2s.reshape(N, -1).astype(np.float32)
    pca = PCA(n_components=min(256, ntr - 1), svd_solver="randomized", random_state=SEED)
    Xtr = pca.fit_transform(X[:ntr]); Xte = pca.transform(X[ntr:])
    k_opt = np.array([kappa_feats(np.stack([unitary_log(Uo[i, k]) for k in range(2)])) for i in range(N)])
    k_rel = np.array([kappa_feats(dk[i]) for i in range(N)])
    Hf = np.array([herm_feats(h_encoding(Uo[i])) for i in range(N)])
    zsel = np.triu(mask)
    Zf = np.array([Zo[i][:, zsel].ravel() for i in range(N)])
    dZf = np.array([dZ[i][:, zsel].ravel() for i in range(N)])
    Zif = np.array([Zi_m[i][:, zsel].ravel() for i in range(N)])
    th_idx = rng.choice(t2hat[0].size, 500, replace=False)
    t2hf = t2hat.reshape(N, -1)[:, th_idx]
    targets = [
        ("kappa_opt canonical (log U_opt)", k_opt),
        ("dkappa (log U_init^dag U_opt)", k_rel),
        ("H(U_opt) phase-free entries", Hf),
        ("Z_opt (allowed entries)", Zf),
        ("dZ = Z_opt - Z_init", dZf),
        ("t2hat_opt (500 entries)", t2hf),
        ("Z_init masked (exact-DF control)", Zif),
        ("lam0*v0 (invariant baseline)", lamv),
    ]
    preds = {}
    print(f"  {'target':36s} {'R2 pooled':>10s} {'R2 mol-specific':>16s}")
    for name, Y in targets:
        (r2p, r2c), a, pred = fit_ridge(Xtr, Y[:ntr], Xte, Y[ntr:])
        preds[name] = pred
        print(f"  {name:36s} {r2p:+10.3f} {r2c:+16.3f}   (alpha={a}, dim={Y.shape[1]})")

    # gauge-aware check of the H-route: rebuild U from predicted H by eigh and
    # measure phase-aligned distance to the true U_opt (relative), vs baselines
    n = norb
    d_pred, d_mean, d_init = [], [], []
    Hmean = h_encoding(Uo[:ntr]).mean(0)
    for t, i in enumerate(range(ntr, N)):
        hp = preds["H(U_opt) phase-free entries"][t]
        # unpack features -> Hermitian (2,n,n)
        H = np.zeros((2, n, n), dtype=complex)
        iu0 = np.triu_indices(n, 0); iu1 = np.triu_indices(n, 1)
        nre = 2 * len(iu0[0])
        re = hp[:nre].reshape(2, -1); im = hp[nre:].reshape(2, -1)
        for k in range(2):
            H[k][iu0] = re[k]; H[k] = H[k] + H[k].T - np.diag(np.diag(H[k]))
            A = np.zeros((n, n), dtype=complex); A[iu1] = 1j * im[k]
            H[k] = H[k] + A + A.conj().T
        Up = np.stack([np.linalg.eigh(H[k])[1] for k in range(2)])
        Um = np.stack([np.linalg.eigh(Hmean[k])[1] for k in range(2)])
        from .compressed_canonical import phase_dist
        d_pred.append(np.sqrt(sum(phase_dist(Uo[i, k], Up[k]) ** 2 for k in range(2))) / np.sqrt(2 * n))
        d_mean.append(np.sqrt(sum(phase_dist(Uo[i, k], Um[k]) ** 2 for k in range(2))) / np.sqrt(2 * n))
        d_init.append(np.sqrt(sum(phase_dist(Uo[i, k], Ui[i, k]) ** 2 for k in range(2))) / np.sqrt(2 * n))
    print(f"\n  U from predicted H (eigh) -> phase-aligned rel. distance to U_opt:"
          f" median {np.median(d_pred):.3f}   [train-mean H: {np.median(d_mean):.3f};"
          f"  canonical exact init: {np.median(d_init):.3f}]")

    # 1-NN label transfer: raw vs phase-orbit distance
    d2 = pair_sqdist(X)[ntr:, :ntr]
    nn = np.argmin(d2, axis=1)
    raw_nn, orb_nn, raw_rnd, orb_rnd = [], [], [], []
    for t, i in enumerate(range(ntr, N)):
        for src, ro, oo in [(int(nn[t]), raw_nn, orb_nn), (int(rng.integers(0, ntr)), raw_rnd, orb_rnd)]:
            ro.append(np.linalg.norm(Uo[i] - Uo[src]))
            oo.append(np.sqrt(sum(phase_dist(Uo[i, k], Uo[src, k]) ** 2 for k in range(2))))
    print(f"  1-NN transfer of U_opt: RAW   1-NN {np.mean(raw_nn):.3f} vs random {np.mean(raw_rnd):.3f}"
          f" (ratio {np.mean(raw_nn)/np.mean(raw_rnd):.3f})")
    print(f"                          ORBIT 1-NN {np.mean(orb_nn):.3f} vs random {np.mean(orb_rnd):.3f}"
          f" (ratio {np.mean(orb_nn)/np.mean(orb_rnd):.3f})")

    # stratify smoothness/learnability by outer gap (Davis-Kahan)
    print("\n--- learnability of kappa_opt stratified by outer gap (test set) ---")
    gte = gaps[ntr:]
    pk = preds["kappa_opt canonical (log U_opt)"]
    for lo, hi in ((0, 0.05), (0.05, 0.2), (0.2, 10)):
        sel = (gte >= lo) & (gte < hi)
        if sel.sum() > 5:
            r2 = r2_score(k_opt[ntr:][sel].ravel(), pk[sel].ravel())
            print(f"  gap in [{lo},{hi}): n={sel.sum():4d}  R^2={r2:+.3f}")

    if args.save:
        np.savez(args.save, **{k.replace(" ", "_"): v for k, v in cols.items()})
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
