"""Canonical-gauge compressed (optimize=True) double-factorization labels.

ffsim's compressed DF works in two steps: an exact truncated DF *init*
(nested eigh -> LAPACK-arbitrary per-column phases, arbitrary basis inside
degenerate eigenspaces, arbitrary sign of v0) followed by L-BFGS on
kappa = log U and the allowed Z entries.  The reconstruction objective is
exactly invariant under the gauge group, but the *optimizer trajectory in
kappa-space is not*, so the optimized labels inherit -- and scramble -- the
init's arbitrary gauge.  That is the source of the "kappa is unlearnable"
observation for optimize=True labels (docs/eigenvector_learnability_study.md).

This module makes the whole label map t2 -> (U, Z) a deterministic, gauge-
canonical function of t2 by canonicalizing the init *before* optimizing:

  1. sign of the dominant t2 eigenvector v0: largest-|entry| element positive;
  2. degenerate (zero-eigenvalue) subspace of the quadrature Q+/-: replaced by a
     deterministic orthonormal basis (pivoted QR of the gauge-invariant
     projector), removing LAPACK's arbitrary null-space rotation;
  3. per-column phases of U: each column rotated so its largest-modulus entry
     is real and positive ('maxabs'), or so its diagonal entry is real positive
     ('diag');
  4. column order: ascending inner eigenvalue w (eigh order), which is already
     deterministic once (1) fixes the sign.

The compressed optimization then runs from this canonical init with ffsim's own
parameterization and objective (reused from ffsim.linalg.util), including the
diagonal-Coulomb norm regularizer of Lin et al. (arXiv:2511.22476).

Also provides the *relative* (residual) parameterization used by the model:
    U_opt = U_init @ expm(dkappa),   Z_opt = Z_init + dZ
with dkappa = log(U_init^dag U_opt) anti-Hermitian.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.linalg
import scipy.optimize

from .common import outer_spectrum, quadrature, vec_to_onebody


# ----------------------------------------------------------------- helpers

def interaction_pairs(connectivity: str, norb: int):
    """ffsim.variational.util.interaction_pairs_spin_balanced, re-implemented
    (avoids importing ffsim for the label pipeline's structural code)."""
    if connectivity == "all-to-all":
        return None, None
    if connectivity == "square":
        return [(p, p + 1) for p in range(norb - 1)], [(p, p) for p in range(norb)]
    if connectivity == "hex":
        return ([(p, p + 1) for p in range(norb - 1)],
                [(p, p) for p in range(norb) if p % 2 == 0])
    if connectivity == "heavy-hex":
        return ([(p, p + 1) for p in range(norb - 1)],
                [(p, p) for p in range(norb) if p % 4 == 0])
    raise ValueError(connectivity)


def diag_coulomb_indices_for(connectivity: str, norb: int):
    pairs_aa, pairs_ab = interaction_pairs(connectivity, norb)
    if pairs_aa is None and pairs_ab is None:
        return None
    return sorted(set((pairs_aa or []) + (pairs_ab or [])))


def z_mask(connectivity: str, norb: int) -> np.ndarray:
    idx = diag_coulomb_indices_for(connectivity, norb)
    if idx is None:
        return np.ones((norb, norb), dtype=bool)
    m = np.zeros((norb, norb), dtype=bool)
    r, c = zip(*idx)
    m[list(r), list(c)] = True
    m[list(c), list(r)] = True
    return m


def canonical_phases(U: np.ndarray, mode: str = "maxabs") -> np.ndarray:
    """Rotate each column of U by a phase so a reference entry is real positive."""
    n = U.shape[1]
    if mode == "maxabs":
        ref = U[np.argmax(np.abs(U), axis=0), np.arange(n)]
    elif mode == "diag":
        ref = np.diag(U)
    else:
        raise ValueError(mode)
    ph = np.ones(n, dtype=complex)
    nz = np.abs(ref) > 1e-12
    ph[nz] = ref[nz] / np.abs(ref[nz])
    return U * ph.conj()[None, :]


def canonical_eigh(H: np.ndarray, zero_tol: float = 1e-9, phase_mode: str = "maxabs"):
    """eigh with a deterministic basis for the (near-)zero eigenspace and
    canonical column phases.  Returns (w, U) in ascending-w order."""
    w, U = np.linalg.eigh(H)
    scale = max(np.max(np.abs(w)), 1e-300)
    null = np.abs(w) < zero_tol * scale
    k = int(null.sum())
    if k > 1:
        # gauge-invariant projector onto the null space -> pivoted QR gives a
        # deterministic orthonormal basis of its range.
        V = U[:, null]
        P = V @ V.conj().T
        Q, _, _ = scipy.linalg.qr(P, pivoting=True, mode="economic")
        U = U.copy()
        U[:, null] = Q[:, :k]
        w = w.copy()
        w[null] = 0.0
    return w, canonical_phases(U, phase_mode)


@dataclass
class ExactInit:
    lam0: float
    v0: np.ndarray                 # (nocc*nvirt,) unit, sign-canonical
    gap: float                     # (|lam0|-|lam1|)/|lam0|
    Z: np.ndarray                  # (2, norb, norb) rank-1, +/- lam0 w w^T
    U: np.ndarray                  # (2, norb, norb) complex, canonical gauge
    w: np.ndarray                  # (2, norb) inner eigenvalues
    znorm_full: float              # sum_k ||Z_k||_F^2 over the FULL exact DF
    outer_eigs: np.ndarray


def canonical_exact_init(t2: np.ndarray, phase_mode: str = "maxabs",
                         n_reps: int = 2) -> ExactInit:
    """Canonical-gauge version of ffsim's optimize=False n_reps=2 init."""
    assert n_reps == 2, "sufficient statistic (lam0, v0) is for n_reps=2"
    nocc, _, nvirt, _ = t2.shape
    eigs, vecs = outer_spectrum(t2)
    lam0 = float(eigs[0])
    v0 = vecs[:, 0]
    v0 = v0 * np.sign(v0[np.argmax(np.abs(v0))])
    gap = (abs(eigs[0]) - abs(eigs[1])) / abs(eigs[0]) if len(eigs) > 1 else 1.0
    M = vec_to_onebody(v0, nocc, nvirt)
    Zs, Us, ws = [], [], []
    for sign, coeff in ((1, lam0), (-1, -lam0)):
        w, U = canonical_eigh(quadrature(M, sign), phase_mode=phase_mode)
        Zs.append(coeff * np.outer(w, w))
        Us.append(U)
        ws.append(w)
    # ffsim's regularizer anchors sum ||Z||^2 to the FULL exact factorization:
    # each outer eigenpair contributes two reps with ||Z||_F^2 = lam^2 ||w||^4,
    # and ||w||^2 = trace(Q^2) = ||M||_F^2 = 1  ->  2 * sum_k lam_k^2.
    # (ffsim truncates the outer spectrum at tol=1e-8 on cumulative |lam|,
    # a negligible difference.)
    znorm_full = 2.0 * float(np.sum(eigs ** 2))
    return ExactInit(lam0=lam0, v0=v0, gap=float(gap), Z=np.array(Zs),
                     U=np.array(Us), w=np.array(ws), znorm_full=znorm_full,
                     outer_eigs=eigs)


# ----------------------------------------------------- compressed optimizer

def reconstruct_t2(Z: np.ndarray, U: np.ndarray, nocc: int) -> np.ndarray:
    full = 1j * np.einsum("kpq,kap,kip,kbq,kjq->ijab", Z, U, U.conj(), U, U.conj(),
                          optimize=True)
    return full[:nocc, :nocc, nocc:, nocc:].real


def rel_residual(t2: np.ndarray, Z: np.ndarray, U: np.ndarray) -> float:
    nocc = t2.shape[0]
    return float(np.linalg.norm(reconstruct_t2(Z, U, nocc) - t2) / np.linalg.norm(t2))


def compress_from_init(
    t2: np.ndarray,
    Z0: np.ndarray,
    U0: np.ndarray,
    *,
    connectivity: str = "all-to-all",
    regularization: float = 0.0,
    znorm_ref: float | None = None,
    maxiter: int = 500,
    method: str = "L-BFGS-B",
    x64: bool = True,
):
    """ffsim's compressed-DF optimization (same parameterization, objective and
    regularizer) started from an arbitrary (Z0, U0) init.

    Returns dict(Z, U, x, nit, nfev, success, message, fun, time).
    """
    import jax
    import jax.numpy as jnp
    from opt_einsum import contract
    from ffsim.linalg.util import df_tensors_from_params, df_tensors_from_params_jax, \
        df_tensors_to_params

    if x64:
        jax.config.update("jax_enable_x64", True)

    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    n_tensors = Z0.shape[0]
    dci = diag_coulomb_indices_for(connectivity, norb)
    if znorm_ref is None:
        znorm_ref = float(np.sum(np.abs(Z0) ** 2))
    t2_j = jnp.asarray(t2)

    def fun(x):
        Z, U = df_tensors_from_params_jax(x, n_tensors, norb, dci)
        rec = 1j * contract("mpq,map,mip,mbq,mjq->ijab", Z, U, U.conj(), U, U.conj(),
                            optimize="greedy")[:nocc, :nocc, nocc:, nocc:]
        loss = 0.5 * jnp.sum(jnp.abs(rec - t2_j) ** 2)
        if regularization:
            loss += regularization * jnp.abs(jnp.sum(jnp.abs(Z) ** 2) - znorm_ref)
        return loss

    vg = jax.jit(jax.value_and_grad(fun))

    def vg_np(x):
        v, g = vg(jnp.asarray(x))
        return float(v), np.asarray(g, dtype=float)

    # mask the init's Z to the allowed pattern (ffsim does this implicitly via
    # the parameterization: disallowed entries are simply not parameters)
    x0 = df_tensors_to_params(Z0, U0, dci)
    t0 = time.time()
    res = scipy.optimize.minimize(vg_np, x0, method=method, jac=True,
                                  options=dict(maxiter=maxiter))
    Z, U = df_tensors_from_params(res.x, n_tensors, norb, dci)
    return dict(Z=np.asarray(Z, dtype=float), U=np.asarray(U), x=res.x,
                nit=int(res.nit), nfev=int(res.nfev), success=bool(res.success),
                message=str(res.message), fun=float(res.fun),
                time=time.time() - t0)


# ------------------------------------------------ relative parameterization

def unitary_log(U: np.ndarray) -> np.ndarray:
    """Principal log of a unitary via complex Schur (robust)."""
    T, Q = scipy.linalg.schur(U.astype(complex), output="complex")
    return Q @ np.diag(1j * np.angle(np.diag(T))) @ Q.conj().T


def relative_kappa(U_init: np.ndarray, U_opt: np.ndarray) -> np.ndarray:
    """dkappa with U_opt = U_init @ expm(dkappa), per rep."""
    return np.stack([unitary_log(U_init[k].conj().T @ U_opt[k])
                     for k in range(U_init.shape[0])])


def phase_align(U_ref: np.ndarray, U: np.ndarray) -> np.ndarray:
    """Closed-form per-column phase alignment of U onto U_ref (no permutation)."""
    ov = np.einsum("ij,ij->j", U_ref.conj(), U)
    ph = np.ones_like(ov)
    nz = np.abs(ov) > 1e-12
    ph[nz] = ov[nz] / np.abs(ov[nz])
    return U * ph.conj()[None, :]


def phase_dist(U_ref: np.ndarray, U: np.ndarray) -> float:
    return float(np.linalg.norm(U_ref - phase_align(U_ref, U)))
