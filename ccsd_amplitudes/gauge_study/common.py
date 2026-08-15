"""Shared utilities for the U/J learnability study.

Core objects
------------
exact_df_t2   : numpy reimplementation of ffsim's exact (optimize=False)
                double factorization of t2, returning the same (Z, U) targets
                as generate_df_targets.py plus the underlying spectra.
align_unitaries: optimal gauge alignment (per-column phase + column
                permutation) between two unitaries, via the Hungarian method.
gauge_distance : orbit distance between (U, Z) label pairs.

Label math (matches ffsim.linalg._double_factorized_t2_explicit)
----------------------------------------------------------------
T[(i,a),(j,b)] = t2[i,j,a,b]  (real symmetric, nocc*nvirt square)
eigh(T) -> (lam_k, v_k), sorted by |lam| descending, truncated to top-1
M = mat(v_0)  (norb x norb, virtual-row/occupied-column block)
Q_pm = 0.5*(1 -/+ 1j)*(M +/- 1j*M.T)          (two Hermitian quadratures)
eigh(Q_pm) -> (w_pm, U_pm)                     (U = the kappa targets)
Z_pm = +/- lam_0 * outer(w_pm, w_pm)           (rank-1 diagonal Coulomb)
reps stored = [(+), (-)]; Q_- = conj(Q_+) so U_1 ~ conj(U_0), Z_1 = -Z_0.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm
from scipy.optimize import linear_sum_assignment


# ---------------------------------------------------------------- exact DF

def t2_matrix(t2: np.ndarray) -> np.ndarray:
    """Reshape t2 (nocc,nocc,nvirt,nvirt) -> symmetric (nocc*nvirt)^2 matrix."""
    nocc, _, nvirt, _ = t2.shape
    return t2.transpose(0, 2, 1, 3).reshape(nocc * nvirt, nocc * nvirt)


def outer_spectrum(t2: np.ndarray):
    """Eigendecomposition of the t2 matrix, sorted by |eigenvalue| descending."""
    eigs, vecs = np.linalg.eigh(t2_matrix(t2))
    order = np.argsort(np.abs(eigs))[::-1]
    return eigs[order], vecs[:, order]


def vec_to_onebody(v: np.ndarray, nocc: int, nvirt: int) -> np.ndarray:
    """Outer eigenvector -> one-body matrix M with M[nocc+a, i] = v[i*nvirt+a]."""
    norb = nocc + nvirt
    M = np.zeros((norb, norb))
    M[nocc:, :nocc] = v.reshape(nocc, nvirt).T
    return M


def quadrature(M: np.ndarray, sign: int) -> np.ndarray:
    """ffsim's Hermitian quadrature: 0.5*(1 - sign*1j)*(M + sign*1j*M.T)."""
    return 0.5 * (1 - sign * 1j) * (M + sign * 1j * M.T.conj())


def exact_df_t2(t2: np.ndarray, n_reps: int = 2):
    """Reimplementation of ffsim's exact truncated DF (optimize=False).

    Returns dict with Z (n_reps,norb,norb), U (n_reps,norb,norb) complex,
    outer_eigs (all), inner_eigs (n_reps,norb), v0, lam0.
    """
    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    outer_eigs, outer_vecs = outer_spectrum(t2)

    Zs, Us, ws = [], [], []
    k = 0
    while len(Zs) < n_reps:
        lam, v = outer_eigs[k], outer_vecs[:, k]
        M = vec_to_onebody(v, nocc, nvirt)
        for sign, coeff in ((1, lam), (-1, -lam)):
            if len(Zs) == n_reps:
                break
            w, U = np.linalg.eigh(quadrature(M, sign))
            Zs.append(coeff * np.outer(w, w))
            Us.append(U)
            ws.append(w)
        k += 1
    return dict(
        Z=np.array(Zs), U=np.array(Us), inner_eigs=np.array(ws),
        outer_eigs=outer_eigs, v0=outer_vecs[:, 0], lam0=outer_eigs[0],
    )


def invariant_B(t2: np.ndarray) -> np.ndarray:
    """Gauge-invariant content of the n_reps=2 exact-DF target:
    the best rank-1 approximation lam0 * v0 v0^T of the t2 matrix,
    returned in t2 tensor layout (nocc,nocc,nvirt,nvirt)."""
    nocc, _, nvirt, _ = t2.shape
    eigs, vecs = outer_spectrum(t2)
    B = eigs[0] * np.outer(vecs[:, 0], vecs[:, 0])
    return B.reshape(nocc, nvirt, nocc, nvirt).transpose(0, 2, 1, 3)


def reconstruct_t2(Z: np.ndarray, U: np.ndarray, nocc: int) -> np.ndarray:
    """t2-hat from (Z, U) stacks; same einsum as ffsim.linalg.reconstruct_t2."""
    full = 1j * np.einsum(
        "kpq,kap,kip,kbq,kjq->ijab", Z, U, U.conj(), U, U.conj(), optimize=True
    )
    return full[:nocc, :nocc, nocc:, nocc:]


# ------------------------------------------------------------- alignment

def align_unitaries(Ua: np.ndarray, Ub: np.ndarray, permute: bool = True):
    """Optimal gauge alignment of Ub's columns to Ua's.

    Solves  min_{P perm, D diag-phase} ||Ua - Ub P D||_F  exactly:
    permutation by the Hungarian method on -|Ua^dag Ub|, then the phase of
    each matched column in closed form.

    Returns (Ub_aligned, perm, dist). With permute=False only phases are fixed.
    """
    n = Ua.shape[0]
    M = Ua.conj().T @ Ub                     # M[p, q] = <ua_p, ub_q>
    if permute:
        row, col = linear_sum_assignment(-np.abs(M))
        perm = np.empty(n, dtype=int)
        perm[row] = col                      # column perm[p] of Ub matches p of Ua
    else:
        perm = np.arange(n)
    overlaps = M[np.arange(n), perm]
    phases = np.ones(n, dtype=complex)
    nz = np.abs(overlaps) > 1e-12
    phases[nz] = overlaps[nz] / np.abs(overlaps[nz])
    Ub_aligned = Ub[:, perm] * phases[None, :].conj()
    dist = np.linalg.norm(Ua - Ub_aligned)
    return Ub_aligned, perm, dist


def permute_Z(Z: np.ndarray, perm: np.ndarray) -> np.ndarray:
    """Apply the column permutation of U to Z: Z'[p,q] = Z[perm[p], perm[q]]."""
    return Z[np.ix_(perm, perm)]


def gauge_distance(Ua, Za, Ub, Zb, allow_rep_swap: bool = True):
    """Orbit distance between two label stacks (n_reps, norb, norb).

    Aligns each rep of b to the corresponding rep of a (phase + permutation),
    optionally minimizing over swapping b's reps. Returns dict with aligned
    U/Z distances and the raw distances for comparison.
    """
    n_reps = Ua.shape[0]
    orders = [list(range(n_reps))]
    if allow_rep_swap and n_reps == 2:
        orders.append([1, 0])

    best = None
    for order in orders:
        dU2 = dZ2 = 0.0
        for r, rb in enumerate(order):
            _, perm, d = align_unitaries(Ua[r], Ub[rb])
            dU2 += d ** 2
            dZ2 += np.linalg.norm(Za[r] - permute_Z(Zb[rb], perm)) ** 2
        cand = (np.sqrt(dU2), np.sqrt(dZ2))
        if best is None or cand[0] < best[0]:
            best = cand
    raw_U = np.linalg.norm(Ua - Ub)
    raw_Z = np.linalg.norm(Za - Zb)
    return dict(dU_aligned=best[0], dZ_aligned=best[1], dU_raw=raw_U, dZ_raw=raw_Z)


def kappa_from_U(U: np.ndarray) -> np.ndarray:
    """Principal matrix log (works on a single matrix or a stack)."""
    if U.ndim == 2:
        return logm(U)
    return np.stack([logm(u) for u in U])


# ---------------------------------------------------------------- loading

def load_cache(cache_dir="rhf_cache", norb: int | None = None):
    """Load the 42-molecule compressed-DF cache (optimize=True, LUCJ pairs).

    The cache stores kappa and J (n_reps, 2, norb, norb); U is reconstructed
    as expm(kappa) and Z is the alpha-alpha diagonal-Coulomb block J[:, 0].
    """
    rows = []
    for f in sorted(Path(cache_dir).glob("*.npz")):
        s = np.load(f)
        if "kappa_real" not in s:
            continue
        if norb is not None and int(s["norb"]) != norb:
            continue
        kappa = (np.asarray(s["kappa_real"], dtype=np.float64)
                 + 1j * np.asarray(s["kappa_imag"], dtype=np.float64))
        rows.append(dict(
            name=f.stem,
            t2=np.asarray(s["t2"], dtype=np.float64),
            U=np.stack([expm(k) for k in kappa]),
            Z=np.asarray(s["J"], dtype=np.float64)[:, 0],
            kappa=kappa,
            nocc=int(s["nocc"]), nvirt=int(s["nvirt"]), norb=int(s["norb"]),
        ))
    return rows


def load_n29(data_dir="rhf_dataset_n29", targets_dir="rhf_targets",
             nocc: int = 16, nvirt: int = 13, max_mols: int | None = None,
             with_U: bool = True):
    """Load the uniform-size subset (norb=29) with its exact-DF targets.

    U is reconstructed as expm(kappa) from the stored float32 kappas.
    """
    rows = []
    files = sorted(Path(data_dir).glob("*.npz"))
    for f in files:
        d = np.load(f)
        if int(d["nocc"]) != nocc or int(d["nvirt"]) != nvirt:
            continue
        tf = Path(targets_dir) / f.name
        if not tf.exists():
            continue
        t = np.load(tf)
        row = dict(
            name=f.stem,
            t2=np.asarray(d["t2"], dtype=np.float64),
            t1=np.asarray(d["t1"], dtype=np.float64),
            Z=np.asarray(t["Z"], dtype=np.float64),
            kappa=np.asarray(t["kappa_real"], dtype=np.float64)
            + 1j * np.asarray(t["kappa_imag"], dtype=np.float64),
            nocc=nocc, nvirt=nvirt, norb=nocc + nvirt,
        )
        if with_U:
            row["U"] = np.stack([expm(k) for k in row["kappa"]])
        rows.append(row)
        if max_mols is not None and len(rows) >= max_mols:
            break
    return rows


def nearest_decile_ratio(d_in: np.ndarray, d_out: np.ndarray, frac=0.1):
    """The deck's metric: mean output distance among the closest-in-input
    decile of pairs, relative to the mean over all pairs."""
    k = max(5, int(len(d_in) * frac))
    idx = np.argsort(d_in)[:k]
    return d_out[idx].mean() / d_out.mean()
