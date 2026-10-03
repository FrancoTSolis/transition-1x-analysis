"""Variable projection for the masked compressed-DF objective: closed-form optimal Z for a given U.

The complex reconstruction is linear in the free (masked, symmetric) Z entries z_f, f = (k, p<=q):
    rec = sum_f z_f B_f,   B_f = i (W_kpq + W_kqp)  (p != q),   B_f = i W_kpp  (p == q)
with W_kpq[i,j,a,b] = X_kp[a,i] X_kq[b,j] and X_kp[a,i] = U_k[a,p] conj(U_k[i,p]).
Everything reduces to n x n objects (Gram formulation):
    <X_kp, X_k'r> = M_kk'[p,r] = (V_k^H V_k')[p,r] * conj((O_k^H O_k')[p,r])   (V/O = virtual/occupied rows)
    <W_kpq, W_k'rs> = M_kk'[p,r] M_kk'[q,s]
    <W_kpq, t2>     = G_k[p,q] = (conj(Xf_k) T2 conj(Xf_k)^T)[p,q],  T2[(a,i),(b,j)] = t2[i,j,a,b]
ffsim's objective  0.5 ||rec - t2||^2 + lam |z^T D z - c|  (D = 1 diag, 2 off-diag; c = sum ||Z_exact||_F^2)
is on the branch z^T D z >= c the ridge problem  (Re(A^H A) + 2 lam D) z = Re(A^H t2),  solved here in float64.
By the envelope theorem the gradient of min_Z F(U, Z) w.r.t. U equals the partial gradient at Z = Z*(U).
"""
from __future__ import annotations

import torch
from torch import Tensor


def free_pairs(mask: Tensor):
    """Upper-triangular allowed (p, q) pairs of a symmetric 0/1 mask -> (P, Q)."""
    iu = torch.nonzero(torch.triu(mask) > 0)
    return iu[:, 0], iu[:, 1]


def slot_vectors(U: Tensor, nocc: int) -> Tensor:
    """Xf (B, R, n, nv*no): Xf[k, p, (a,i)] = U_k[a,p] conj(U_k[i,p])."""
    B, R, n, _ = U.shape
    occ, vir = U[..., :nocc, :], U[..., nocc:, :]
    X = vir[:, :, :, None, :] * occ.conj()[:, :, None, :, :]          # (B,R,nv,no,n)
    return X.reshape(B, R, -1, n).transpose(-1, -2)                     # (B,R,n,nv*no)


def t2_matrix_ai(t2: Tensor) -> Tensor:
    """T2[(a,i),(b,j)] = t2[i,j,a,b]  (B, nv*no, nv*no)."""
    B, no, _, nv, _ = t2.shape
    return t2.permute(0, 3, 1, 4, 2).reshape(B, nv * no, nv * no)


def gram(U: Tensor, nocc: int) -> Tensor:
    """M[b, k, k', p, r] = <X_kp, X_k'r>  (B, R, R, n, n) complex."""
    O, V = U[..., :nocc, :], U[..., nocc:, :]
    VV = torch.einsum("bkap,blar->bklpr", V.conj(), V)
    OO = torch.einsum("bkip,blir->bklpr", O.conj(), O)
    return VV * OO.conj()


def project(Y: Tensor, U: Tensor, nocc: int) -> Tensor:
    """G[b, k, p, q] = <W_kpq, Y> for a (real or complex) 4-index tensor Y (B,no,no,nv,nv)."""
    Xf = slot_vectors(U, nocc)                                           # (B,R,n,N)
    Ym = t2_matrix_ai(Y.to(Xf.dtype) if not torch.is_complex(Y) else Y)
    return torch.einsum("bkpx,bxy,bkqy->bkpq", Xf.conj(), Ym.to(Xf.dtype), Xf.conj())


def solve_z(U: Tensor, t2: Tensor, mask: Tensor, lam: float, znorm_ref: Tensor | None = None,
            P=None, Q=None, ridge_floor: float = 1e-12, return_info: bool = False):
    """Optimal masked symmetric Z (B, R, n, n) for the given U under ffsim's complex objective."""
    B, R, n, _ = U.shape
    nocc = t2.shape[1]
    if P is None:
        P, Q = free_pairs(mask)
    m = len(P)
    U64 = U.to(torch.complex128)
    M = gram(U64, nocc)                                                  # (B,R,R,n,n)
    G = project(t2.to(torch.float64), U64, nocc)                          # (B,R,n,n)
    offd = (P != Q).to(torch.float64)
    # ordered pair lists per free index: (P,Q) weight 1, (Q,P) weight offd
    a1, b1, a2, b2 = P, Q, Q, P
    # AhA[k f, l g] = sum_{alpha,beta} w_alpha(f) w_beta(g) M_kl[a_alpha(f), a_beta(g)] M_kl[b_alpha(f), b_beta(g)]
    def blk(aa, bb, cc, dd):
        return M[:, :, :, aa][:, :, :, :, cc] * M[:, :, :, bb][:, :, :, :, dd]   # (B,R,R,m,m)
    A11 = blk(a1, b1, a1, b1)
    A12 = blk(a1, b1, a2, b2) * offd[None, None, None, None, :]
    A21 = blk(a2, b2, a1, b1) * offd[None, None, None, :, None]
    A22 = blk(a2, b2, a2, b2) * (offd[:, None] * offd[None, :])[None, None, None]
    AhA = (A11 + A12 + A21 + A22).real                                   # (B,R,R,m,m)
    AhA = AhA.permute(0, 1, 3, 2, 4).reshape(B, R * m, R * m)
    # Re(A^H t2)[k f] = Re(-i * (G_k[a1,b1] + offd * G_k[a2,b2])) = Im(...)
    Aht = (G[:, :, a1, b1] + offd * G[:, :, a2, b2]).imag.reshape(B, R * m)
    dvec = torch.where(P == Q, 1.0, 2.0).to(torch.float64).repeat(R).to(U.device)
    eye_f = ridge_floor * torch.eye(R * m, device=U.device, dtype=torch.float64)

    def z_of(mu):            # mu: (B,) multiplier on z^T D z;  minimizer of ||Az-t||^2/2 + (mu/2) z^T D z
        return torch.linalg.solve(AhA + torch.diag_embed(mu[:, None] * dvec[None, :]) + eye_f, Aht)

    def nrm(z):
        return (z * dvec * z).sum(1)

    if lam and znorm_ref is not None:
        # exact minimizer of 0.5||Az-t||^2 + lam |z^T D z - c|: stationarity (AhA + mu D) z = Aht with
        # mu = 2 lam (norm > c), mu = -2 lam (norm < c) or |mu| <= 2 lam with norm = c (kink).
        # norm(z(mu)) decreases monotonically in mu, so: mu=+2lam if norm(z(+2lam)) >= c; mu=-2lam if
        # norm(z(-2lam)) <= c; otherwise bisect mu in (-2lam, 2lam) for norm = c.
        c = znorm_ref.to(torch.float64)
        lo = torch.full((B,), -2.0 * lam, device=U.device, dtype=torch.float64)
        hi = torch.full((B,), 2.0 * lam, device=U.device, dtype=torch.float64)
        with torch.no_grad():
            n_hi, n_lo = nrm(z_of(hi)), nrm(z_of(lo))
            mu = torch.where(n_hi >= c, hi, torch.where(n_lo <= c, lo, torch.zeros_like(hi)))
            kink = (n_hi < c) & (n_lo > c)
            if kink.any():
                a, b = lo.clone(), hi.clone()
                for _ in range(40):
                    mid = 0.5 * (a + b)
                    nm = nrm(z_of(mid))
                    up = nm > c                    # norm too large -> increase mu
                    a = torch.where(up, mid, a)
                    b = torch.where(up, b, mid)
                mu = torch.where(kink, 0.5 * (a + b), mu)
        z = z_of(mu)                               # differentiable in U at fixed mu
    else:
        z = z_of(torch.full((B,), 2.0 * lam, device=U.device, dtype=torch.float64))
    zz = z.reshape(B, R, m)
    Z = torch.zeros(B, R, n, n, device=U.device, dtype=torch.float64)
    Z = Z.index_put((torch.arange(B, device=U.device)[:, None, None], torch.arange(R, device=U.device)[None, :, None],
                     P[None, None, :], Q[None, None, :]), zz)
    Z = Z + Z.transpose(-1, -2) - torch.diag_embed(torch.diagonal(Z, dim1=-2, dim2=-1))
    Z = Z.to(t2.dtype)
    if return_info:
        info = {}
        if znorm_ref is not None and lam:
            info["branch_plus_frac"] = float((mu > 1.999 * lam).float().mean())
            info["branch_minus_frac"] = float((mu < -1.999 * lam).float().mean())
            info["kink_frac"] = float(kink.float().mean())
        return Z, info
    return Z
