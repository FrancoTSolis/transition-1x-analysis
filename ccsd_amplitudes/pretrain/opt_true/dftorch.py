"""Batched torch implementation of ffsim's compressed double-factorized t2 objective.

Shared toolkit for the optimize=True pretraining experiments (Oct 2026). All n29
molecules share one shape (nocc=16, nvirt=13, norb=29), so the whole group can be
reconstructed / optimized as one batch on a GPU.

Conventions (identical to ffsim.linalg.double_factorized_t2, optimize=True):
    rec[i,j,a,b]   = i * sum_k sum_pq Z_k[p,q] U_k[a,p] U_k*[i,p] U_k[b,q] U_k*[j,q]   (complex)
    t2hat          = Re(rec)   (what ffsim's reconstruct_t2(...).real and our "resid" report)
    loss           = 0.5 * ||rec - t2||^2 + lam * | sum_k ||Z_k||_F^2 - znorm_ref |
NOTE: ffsim's loss uses the COMPLEX reconstruction (|rec - t2|^2 includes Im(rec)^2). The exact-DF
init is a conjugate-pair-symmetric saddle (rep1 = conj rep0, Z1 = -Z0, rec real); a real-part-only
loss keeps the flow on that symmetric manifold (best 0.415), the complex loss leaves it (0.32-0.34).
with Z_k real symmetric, nonzero only on the connectivity mask (square: (p,p), (p,p+1)),
and U_k unitary. Gauge: per-column phases of U_k leave t2hat and the LUCJ operator
invariant (Z unchanged).

Parameterization for optimization / prediction:
    U = U_base @ expm(K),  K = A + i S,  A real antisymmetric, S real symmetric
    Z = (Z_base + D) * mask, D real symmetric
"""
from __future__ import annotations

import math

import torch
from torch import Tensor


# ---------------------------------------------------------------- structure

def square_mask(norb: int, device=None, dtype=torch.float32) -> Tensor:
    """Connectivity mask of the 'square' LUCJ pattern: (p,p) and (p,p+/-1)."""
    m = torch.eye(norb, device=device, dtype=dtype)
    idx = torch.arange(norb - 1, device=device)
    m[idx, idx + 1] = 1
    m[idx + 1, idx] = 1
    return m


def full_mask(norb: int, device=None, dtype=torch.float32) -> Tensor:
    return torch.ones(norb, norb, device=device, dtype=dtype)


def antihermitian(A: Tensor, S: Tensor) -> Tensor:
    """K = antisym(A) + i * sym(S) from unconstrained real tensors (..., n, n)."""
    Aa = 0.5 * (A - A.transpose(-1, -2))
    Ss = 0.5 * (S + S.transpose(-1, -2))
    return torch.complex(Aa, Ss)


def sym(D: Tensor) -> Tensor:
    return 0.5 * (D + D.transpose(-1, -2))


def unitary_from(U_base: Tensor, A: Tensor, S: Tensor) -> Tensor:
    """U_base @ expm(antihermitian(A, S)); shapes (..., n, n)."""
    K = antihermitian(A, S).to(U_base.dtype)
    return U_base @ torch.linalg.matrix_exp(K)


# ----------------------------------------------------------- reconstruction

def reconstruct(Z: Tensor, U: Tensor, nocc: int) -> Tensor:
    """Batched t2hat. Z: (B, R, n, n) real, U: (B, R, n, n) complex -> (B, no, no, nv, nv) real."""
    B, R, n, _ = U.shape
    nv = n - nocc
    occ = U[..., :nocc, :]                                   # (B,R,no,n)
    vir = U[..., nocc:, :]                                   # (B,R,nv,n)
    Apair = vir[:, :, :, None, :] * occ.conj()[:, :, None, :, :]   # (B,R,nv,no,n)  [a,i,p]
    Apair = Apair.reshape(B, R, nv * nocc, n)
    Bm = Apair @ Z.to(U.dtype)                               # (B,R,nv*no,n)
    T = torch.einsum("brxq,bryq->bxy", Bm, Apair)            # (B, nv*no, nv*no)
    t2h = -T.imag.reshape(B, nv, nocc, nv, nocc).permute(0, 2, 4, 1, 3)
    return t2h


def reconstruct_complex(Z: Tensor, U: Tensor, nocc: int) -> Tensor:
    """Batched complex rec = i*T (ffsim's `full` before taking .real), (B, no, no, nv, nv)."""
    B, R, n, _ = U.shape
    nv = n - nocc
    occ = U[..., :nocc, :]
    vir = U[..., nocc:, :]
    Apair = (vir[:, :, :, None, :] * occ.conj()[:, :, None, :, :]).reshape(B, R, nv * nocc, n)
    T = torch.einsum("brxq,bryq->bxy", Apair @ Z.to(U.dtype), Apair)
    return (1j * T).reshape(B, nv, nocc, nv, nocc).permute(0, 2, 4, 1, 3)


def rel_residual_imag(Z: Tensor, U: Tensor, t2: Tensor) -> Tensor:
    """||Im rec|| / ||t2|| per molecule (0 for conjugate-paired solutions)."""
    rec = reconstruct_complex(Z, U, t2.shape[1])
    return rec.imag.flatten(1).norm(dim=1) / t2.flatten(1).norm(dim=1)


def rel_residual(Z: Tensor, U: Tensor, t2: Tensor) -> Tensor:
    """||t2hat - t2|| / ||t2|| per molecule (B,)."""
    nocc = t2.shape[1]
    d = reconstruct(Z, U, nocc) - t2
    return d.flatten(1).norm(dim=1) / t2.flatten(1).norm(dim=1)


def objective(Z: Tensor, U: Tensor, t2: Tensor, lam: float = 0.0,
              znorm_ref: Tensor | None = None, reduce: str = "none",
              complex_obj: bool = True) -> Tensor:
    """ffsim's compressed-DF loss per molecule: 0.5||rec-t2||^2 + lam*|sum||Z||^2 - ref|.
    complex_obj=False drops the Im(rec)^2 part (NOT ffsim's objective; see module docstring)."""
    nocc = t2.shape[1]
    if complex_obj:
        d = reconstruct_complex(Z, U, nocc) - t2
        loss = 0.5 * d.abs().pow(2).flatten(1).sum(1)
    else:
        d = reconstruct(Z, U, nocc) - t2
        loss = 0.5 * d.flatten(1).pow(2).sum(1)
    if lam:
        zn = Z.flatten(1).pow(2).sum(1)
        loss = loss + lam * (zn - znorm_ref).abs()
    if reduce == "mean":
        return loss.mean()
    if reduce == "sum":
        return loss.sum()
    return loss


def normalized_objective(Z, U, t2, lam=0.0, znorm_ref=None, complex_obj: bool = True) -> Tensor:
    """Objective divided by 0.5||t2||^2 (scale-free across molecules): |r|^2/||t2||^2 + 2 lam |.|/||t2||^2."""
    t2n2 = t2.flatten(1).pow(2).sum(1)
    return objective(Z, U, t2, lam, znorm_ref, complex_obj=complex_obj) / (0.5 * t2n2)


# ------------------------------------------------------------------ gauge

def column_overlaps(Ua: Tensor, Ub: Tensor) -> Tensor:
    """|<ua_p, ub_p>| per column, (..., n)."""
    return (Ua.conj() * Ub).sum(-2).abs()


def phase_dist(Ua: Tensor, Ub: Tensor) -> Tensor:
    """min over per-column phases of ||Ua - Ub D||_F, summed over reps -> (B,).
    ||u - v e^{it}||^2 minimized = 2 - 2|<u,v>| for unit columns."""
    ov = column_overlaps(Ua, Ub)
    return (2.0 - 2.0 * ov).clamp_min(0).flatten(1).sum(1).sqrt()


def phase_invariant_u_loss(U_pred: Tensor, U_tgt: Tensor) -> Tensor:
    """Gauge-invariant (per-column phase) U loss in [0, 1]: mean_p (1 - |<u_p, v_p>|^2), per molecule."""
    ov2 = column_overlaps(U_pred, U_tgt).pow(2)
    return (1.0 - ov2).flatten(1).mean(1)


# ------------------------------------------------------------ refinement

@torch.no_grad()
def _clone(x):
    return x.detach().clone()


def refine(t2: Tensor, U_base: Tensor, Z_base: Tensor, mask: Tensor, *, lam: float = 0.0,
           znorm_ref: Tensor | None = None, steps: int = 1000, lr: float = 0.02,
           prox_mu: float = 0.0, A0: Tensor | None = None, S0: Tensor | None = None,
           D0: Tensor | None = None, log_every: int = 0, chunk: int | None = None,
           betas=(0.9, 0.99), lr_decay: float = 0.0, return_params: bool = False,
           complex_obj: bool = True):
    """Batched Adam on (A, S, D) from (U_base, Z_base) [+ optional starting params].

    Minimizes the normalized ffsim objective per molecule (sum over the batch, so each
    molecule's parameters get their own Adam statistics). prox_mu adds
    prox_mu * ||(A,S,D) - (A0,S0,D0)||^2 (a trust region around the start).
    Returns (U, Z, history) with history = list of (step, median resid).
    """
    B, R, n, _ = U_base.shape
    dev = U_base.device
    rdt = torch.float32 if U_base.dtype == torch.complex64 else torch.float64
    A = (A0.clone() if A0 is not None else torch.zeros(B, R, n, n, device=dev, dtype=rdt)).requires_grad_(True)
    S = (S0.clone() if S0 is not None else torch.zeros(B, R, n, n, device=dev, dtype=rdt)).requires_grad_(True)
    D = (D0.clone() if D0 is not None else torch.zeros(B, R, n, n, device=dev, dtype=rdt)).requires_grad_(True)
    cA, cS, cD = _clone(A), _clone(S), _clone(D)
    opt = torch.optim.Adam([A, S, D], lr=lr, betas=betas)
    sched = (torch.optim.lr_scheduler.LambdaLR(opt, lambda s: 1.0 / (1.0 + lr_decay * s))
             if lr_decay else None)
    hist = []
    t2n2 = t2.flatten(1).pow(2).sum(1)
    for it in range(steps + 1):
        U = unitary_from(U_base, A, S)
        Z = (Z_base + sym(D)) * mask
        loss_m = objective(Z, U, t2, lam, znorm_ref, complex_obj=complex_obj) / (0.5 * t2n2)
        if log_every and (it % log_every == 0 or it == steps):
            with torch.no_grad():
                r = rel_residual(Z, U, t2)
                hist.append((it, float(r.median()), float(loss_m.median())))
        if it == steps:
            break
        loss = loss_m.sum()
        if prox_mu:
            loss = loss + prox_mu * ((A - cA).pow(2).flatten(1).sum(1) + (S - cS).pow(2).flatten(1).sum(1)
                                     + (D - cD).pow(2).flatten(1).sum(1)).sum()
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        if sched is not None:
            sched.step()
    with torch.no_grad():
        U = unitary_from(U_base, A, S)
        Z = (Z_base + sym(D)) * mask
    if return_params:
        return U.detach(), Z.detach(), hist, (A.detach(), S.detach(), D.detach())
    return U.detach(), Z.detach(), hist


def refine_chunked(t2, U_base, Z_base, mask, chunk: int = 256, **kw):
    """refine() over chunks of the batch (bounded memory); concatenates results."""
    Us, Zs, hists = [], [], []
    for s in range(0, t2.shape[0], chunk):
        sl = slice(s, s + chunk)
        extra = {k: (v[sl] if isinstance(v, Tensor) and v.dim() > 0 and v.shape[0] == t2.shape[0] else v)
                 for k, v in kw.items()}
        U, Z, h = refine(t2[sl], U_base[sl], Z_base[sl], mask, **extra)
        Us.append(U)
        Zs.append(Z)
        hists.append(h)
    return torch.cat(Us), torch.cat(Zs), hists


def unitary_log(U: Tensor) -> Tensor:
    """Principal log of (batched) unitaries via eigendecomposition of the Hermitian generator.
    For a unitary, U = V diag(e^{i th}) V^dag; log = V diag(i th) V^dag (eigh of the
    Hermitian part is not used: use the Schur-free route through torch.linalg.eig)."""
    w, V = torch.linalg.eig(U)
    th = torch.angle(w)
    Vinv = torch.linalg.inv(V)
    return V @ torch.diag_embed(1j * th.to(V.dtype)) @ Vinv


__all__ = [
    "square_mask", "full_mask", "antihermitian", "sym", "unitary_from", "reconstruct",
    "reconstruct_complex", "rel_residual_imag", "rel_residual", "objective", "normalized_objective", "column_overlaps", "phase_dist",
    "phase_invariant_u_loss", "refine", "refine_chunked", "unitary_log",
]
