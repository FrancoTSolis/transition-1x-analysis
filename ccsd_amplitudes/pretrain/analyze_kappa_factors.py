#!/usr/bin/env python3
"""What actually determines kappa? Empirical battery on the RHF cache.

We know t2-slices don't predict kappa (gauge). Here we test whether kappa is
predictable from / correlated with other quantities:
  - J / Z (diagonal Coulomb)         -- same DF, but eigenvalues vs eigenvectors
  - orbital energies & energy gaps   -- the Fock spectrum that sets the basis
  - t2 (summary)                     -- baseline
both ACROSS molecules (regression R^2) and WITHIN molecules (elementwise corr).
"""
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score


def load(norb_only=True):
    S = [np.load(f) for f in sorted(Path("rhf_cache").glob("*.npz"))]
    if norb_only:
        norb = Counter(int(s["norb"]) for s in S).most_common(1)[0][0]
        S = [s for s in S if int(s["norb"]) == norb]
    return S, int(S[0]["norb"]), int(S[0]["nocc"]), int(S[0]["nvirt"])


def regress(X, y, n_tr, alpha=10.0):
    sc = StandardScaler().fit(X[:n_tr])
    Xs = sc.transform(X)
    r = Ridge(alpha).fit(Xs[:n_tr], y[:n_tr])
    return r2_score(y[n_tr:].ravel(), r.predict(Xs[n_tr:]).ravel())


def main():
    S, norb, nocc, nvirt = load()
    n = len(S); n_tr = int(0.8 * n)
    triu = np.triu_indices(norb, 1)
    print(f"{n} molecules, norb={norb} (nocc={nocc}, nvirt={nvirt})\n")

    # ---- targets ----
    KR = np.array([s["kappa_real"][0] for s in S])          # (n, norb, norb)
    y_kappa = KR[:, triu[0], triu[1]]                        # (n, n_pairs)

    # ---- feature sets ----
    Zfull = np.array([s["J"][0, 0] for s in S])             # (n, norb, norb)
    feat_Z = Zfull[:, triu[0], triu[1]]
    feat_eps = np.array([s["orb_energies"] for s in S])     # (n, norb)
    # energy-gap matrix (eps_p - eps_q)
    feat_egap = np.array([(s["orb_energies"][:, None] - s["orb_energies"][None, :])[triu]
                          for s in S])
    feat_t2 = np.array([np.concatenate([s["t2"].mean(axis=(2, 3)).ravel(),
                                        s["t2"].mean(axis=(0, 1)).ravel(),
                                        [s["t2"].std(), np.abs(s["t2"]).max()]]) for s in S])

    print("=== ACROSS molecules: predict kappa_re (upper-tri) via Ridge ===")
    print(f"{'features':<28}{'R^2':>8}")
    for name, X in [("t2 (summary)", feat_t2), ("J / Z", feat_Z),
                    ("orbital energies", feat_eps), ("energy gaps eps_p-eps_q", feat_egap),
                    ("J + energies + gaps", np.concatenate([feat_Z, feat_eps, feat_egap], 1)),
                    ("everything + t2", np.concatenate([feat_t2, feat_Z, feat_eps, feat_egap], 1))]:
        print(f"{name:<28}{regress(X, y_kappa, n_tr):>8.3f}")

    # ---- WITHIN molecule: elementwise correlation over all (p,q) entries ----
    print("\n=== WITHIN molecule: corr(kappa_re[p,q], X[p,q]) averaged over molecules ===")
    def within_corr(getX):
        cs = []
        for s in S:
            kr = s["kappa_real"][0]
            X = getX(s)
            a, b = kr.ravel(), X.ravel()
            if a.std() > 1e-9 and b.std() > 1e-9:
                cs.append(np.corrcoef(a, b)[0, 1])
        return np.mean(cs), np.std(cs)

    eps_gap = lambda s: s["orb_energies"][:, None] - s["orb_energies"][None, :]
    for name, getX in [
        ("J / Z[p,q]", lambda s: s["J"][0, 0]),
        ("|J / Z[p,q]|", lambda s: np.abs(s["J"][0, 0])),
        ("energy gap (eps_p-eps_q)", eps_gap),
        ("1/energy gap", lambda s: 1.0 / (eps_gap(s) + np.sign(eps_gap(s) + 1e-9) * 1e-3)),
    ]:
        m, sd = within_corr(getX)
        print(f"  kappa_re vs {name:<26} corr = {m:+.3f} ± {sd:.3f}")
    # |kappa| vs |J| and vs |gap| (magnitudes)
    print("\n  (magnitudes) corr(|kappa_re[p,q]|, .):")
    for name, getX in [("|J/Z|", lambda s: np.abs(s["J"][0, 0])),
                       ("1/|gap|", lambda s: 1.0 / (np.abs(eps_gap(s)) + 1e-3))]:
        cs = []
        for s in S:
            a = np.abs(s["kappa_real"][0]).ravel(); b = getX(s).ravel()
            if a.std() > 1e-9 and b.std() > 1e-9:
                cs.append(np.corrcoef(a, b)[0, 1])
        print(f"    |kappa| vs {name:<10} corr = {np.mean(cs):+.3f} ± {np.std(cs):.3f}")

    # ---- per-molecule: does ||kappa|| track spectral crowding? ----
    print("\n=== PER molecule: does ||kappa|| track near-degeneracy? ===")
    kn = np.array([np.linalg.norm(s["kappa_real"][0]) for s in S])
    # min gap between consecutive sorted orbital energies (spectral crowding)
    mingap = np.array([np.min(np.diff(np.sort(s["orb_energies"]))) for s in S])
    # min gap between sorted Z eigenvalues (DF spectral crowding)
    zmingap = np.array([np.min(np.abs(np.diff(np.sort(np.linalg.eigvalsh(s["J"][0, 0])))))
                        for s in S])
    print(f"  corr(||kappa||, 1/min orbital-energy gap) = "
          f"{np.corrcoef(kn, 1/(mingap+1e-6))[0,1]:+.3f}")
    print(f"  corr(||kappa||, 1/min Z-eigenvalue gap)    = "
          f"{np.corrcoef(kn, 1/(zmingap+1e-9))[0,1]:+.3f}")


if __name__ == "__main__":
    main()
