"""Discrete symmetry group of the masked compressed-DF objective (n_reps = 2).

Beyond per-column phases of each U_k, the objective
    0.5 ||rec - t2||^2 + lam |sum ||Z||^2 - ref|,   rec = i sum_k sum_pq Z_k U_k[a,p] U_k*[i,p] U_k[b,q] U_k*[j,q]
(t2 real) is exactly invariant under
    swap   : (U0,Z0,U1,Z1) -> (U1,Z1,U0,Z0)                 (sum over k)
    conj   : U_k -> conj(U_k), Z_k -> -Z_k                    (rec -> conj(rec); t2 real)
    rev    : U_k -> U_k P, Z_k -> P Z_k P, P = index reversal (maps the square mask to itself)
    gamma  : U_k -> Gamma U_k, Gamma = diag(+1_occ, -1_virt)  (every rec term carries two virtual rows)
The exact-DF init is a fixed point of (conj o swap) (U1 = conj U0, Z1 = -Z0, rec real), i.e. a
symmetric saddle: the compressed optimizer breaks the symmetry spontaneously (driven by round-off),
so each label is a random representative of its orbit under this group. Supervised targets must be
aligned over the group (and phases) before regression, or the loss must be the orbit minimum.

Note: rev and gamma are symmetries of the t2 reconstruction (first-order DF), not of the LUCJ
energy in general; swap reorders the two LUCJ layers (the product is not commutative). Energies of
orbit partners can therefore differ; the alignment is for learning the DF target.
"""
from __future__ import annotations

import itertools

import torch
from torch import Tensor

NOCC = 16


def apply(g: tuple[int, int, int, int], U: Tensor, Z: Tensor, nocc: int = NOCC):
    """Apply group element g = (swap, conj, rev, gamma) in {0,1}^4 to batched (U, Z) (B, 2, n, n)."""
    sw, cj, rv, gm = g
    if sw:
        U = U.flip(1)
        Z = Z.flip(1)
    if cj:
        U = U.conj()
        Z = -Z
    if rv:
        U = U.flip(-1)
        Z = Z.flip(-1).flip(-2)
    if gm:
        s = torch.ones(U.shape[-2], device=U.device, dtype=U.real.dtype)
        s[nocc:] = -1
        U = U * s[:, None].to(U.dtype)
    return U, Z


def elements(swap=True, conj=True, rev=True, gamma=True):
    choices = [(0, 1) if f else (0,) for f in (swap, conj, rev, gamma)]
    return list(itertools.product(*choices))


def phase_sqdist(Ua: Tensor, Ub: Tensor) -> Tensor:
    """sum over reps/columns of min-phase squared distance, (B,) for broadcastable (B,2,n,n)."""
    ov = (Ua.conj() * Ub).sum(-2).abs()
    return (2.0 - 2.0 * ov).clamp_min(0).flatten(1).sum(1)


def orbit_align(U_ref: Tensor, Z_ref: Tensor, U: Tensor, Z: Tensor, group=None, z_weight: float = 0.0):
    """For each molecule pick g minimizing phase_sqdist(U_ref, gU) (+ z_weight ||Z_ref - gZ||^2),
    then fix column phases of gU to U_ref. Returns (U_aligned, Z_aligned, best_g_index, sqdist)."""
    group = group or elements()
    best = None
    for gi, g in enumerate(group):
        Ug, Zg = apply(g, U, Z)
        dist = phase_sqdist(U_ref, Ug)
        if z_weight:
            dist = dist + z_weight * (Z_ref - Zg).flatten(1).pow(2).sum(1)
        if best is None:
            best = (dist, torch.full_like(dist, gi, dtype=torch.long))
        else:
            better = dist < best[0]
            best = (torch.where(better, dist, best[0]), torch.where(better, torch.full_like(best[1], gi), best[1]))
    dist, gidx = best
    Ua = torch.empty_like(U)
    Za = torch.empty_like(Z)
    for gi, g in enumerate(group):
        sel = gidx == gi
        if sel.any():
            Ug, Zg = apply(g, U[sel], Z[sel])
            ov = (U_ref[sel].conj() * Ug).sum(-2)                     # (b,2,n)
            ph = torch.where(ov.abs() > 1e-12, ov / ov.abs().clamp_min(1e-30), torch.ones_like(ov))
            Ua[sel] = Ug * ph.conj()[:, :, None, :]
            Za[sel] = Zg
    return Ua, Za, gidx, dist


def pairwise_orbit_sqdist(U: Tensor, group=None, chunk: int = 64) -> Tensor:
    """(N, N) matrix of min over group (and phases) of squared U distance."""
    group = group or elements()
    N = U.shape[0]
    out = torch.empty(N, N, device=U.device, dtype=U.real.dtype)
    Zdummy = torch.zeros(U.shape, device=U.device, dtype=U.real.dtype)
    variants = [apply(g, U, Zdummy)[0] for g in group]                            # each (N,2,n,n)
    for s in range(0, N, chunk):
        a = U[s:s + chunk]                                                        # (c,2,n,n)
        best = None
        for V in variants:
            # |<a_p, V_p>| for all pairs: (c, N, 2, n)
            ov = (a.conj()[:, None] * V[None]).sum(-2).abs()
            dist = (2.0 - 2.0 * ov).clamp_min(0).flatten(2).sum(2)
            best = dist if best is None else torch.minimum(best, dist)
        out[s:s + chunk] = best
    return out
