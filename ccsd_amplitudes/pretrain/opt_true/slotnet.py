"""Slot-frame learned optimizer for the masked compressed-DF problem (frame-consistent pretraining model).

Diagnosis behind this model: the old heads emit the generator dK in the frame of the init's columns (slots),
U = U0 expm(dK), but read it from tokens indexed by MO pairs; U0's slot basis jumps discontinuously between
molecules (near-degenerate inner eigenvalues), so no continuous MO-frame network can learn the compensation.
Here the tokens ARE slot pairs of the current iterate U_t, and every input feature is expressed in that same
frame, so inputs and outputs jump together.

State: U_t (B, 2, n, n) complex; Z_t = Z*(U_t) by exact variable projection (pretrain.opt_true.varpro).
Features per slot pair (p, q), both reps as channels (all column-phase invariant except the gradient, which is
covariant; column phases are canonicalized every recycle):
    G_k[p,q]  = <W_kpq, t2>/||t2||             (Re, Im)   projection of the target on slot pair (p,q)
    R_k[p,q]  = <W_kpq, t2 - rec>/||t2||       (Re, Im)   residual projected the same way
    g_k[p,q]  = dF/d(dK_k)[p,q] at dK = 0      (Re, Im)   per-molecule RMS-normalized (Danskin, Z* fixed)
    Z*_k[p,q], mask, |p-q| one-hot (0..3+)
    M_kk[p,q] (|.|, Re, Im) and M_01[p,q] (Re, Im)   Gram overlaps <X_kp, X_k'q>
    node: occupied weight o_kp, o_kq; Z*_kpp, Z*_kqq; sinusoidal slot position of p and q
Trunk: Linear -> d, + LayerNorm(previous recycle's pair state), L x TransformerLayer (pretrain.model).
Heads: per rep, dK_k = kscale * tanh(f(z_pq) - f(z_qp)) + i kscale * tanh(h(z_pq) + h(z_qp)), zero-init.
Update: U_{t+1} = canon_phases(U_t expm(dK)).   T = 1 is the frame-aware ONE-SHOT pretraining model.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
from torch import Tensor

from pretrain.model import ModelConfig, TransformerLayer
from pretrain.opt_true import dftorch as D
from pretrain.opt_true import varpro as V

N_FEAT_PAIR = 2 * (2 + 2 + 2 + 1 + 3) + 2 + 1 + 4      # per-rep pair feats x2 reps, cross Gram, mask, |p-q|
N_POS = 8


def canon_phases(U: Tensor) -> Tensor:
    """Rotate each column so its largest-modulus entry is real positive (exact gauge transformation)."""
    idx = U.abs().argmax(dim=-2, keepdim=True)                       # (B,R,1,n)
    ref = torch.gather(U, -2, idx)
    ph = ref / ref.abs().clamp_min(1e-30)
    return U * ph.conj()


def sinusoid(n: int, dim: int, device) -> Tensor:
    pos = torch.arange(n, device=device, dtype=torch.float32)[:, None]
    k = torch.arange(dim // 2, device=device, dtype=torch.float32)[None, :]
    ang = pos / (10.0 ** (k / max(1, dim // 2 - 1)))
    return torch.cat([torch.sin(ang), torch.cos(ang)], -1)            # (n, dim)


@torch.no_grad()
def features(U: Tensor, t2: Tensor, mask: Tensor, lam: float, zref: Tensor):
    """Slot-pair features (B, n, n, F) for state U, plus Z* and the per-molecule normalized objective."""
    U = U.detach()
    B, R, n, _ = U.shape
    nocc = t2.shape[1]
    Z = V.solve_z(U, t2, mask, lam, zref)
    tn = t2.flatten(1).norm(dim=1).to(torch.float64)
    U64 = U.to(torch.complex128)
    rec = D.reconstruct_complex(Z.to(torch.float64), U64, nocc)
    Rres = t2.to(torch.float64) - rec
    G = V.project(t2.to(torch.float64), U64, nocc) / tn[:, None, None, None]
    GR = V.project(Rres, U64, nocc) / tn[:, None, None, None]
    M = V.gram(U64, nocc)                                              # (B,R,R,n,n)
    # gradient wrt right generator at 0 with Z fixed (Danskin): d/dK of F(U expm(K), Z*)
    with torch.enable_grad():
        A = torch.zeros(B, R, n, n, device=U.device, dtype=torch.float32, requires_grad=True)
        S = torch.zeros_like(A, requires_grad=True)
        Uk = D.unitary_from(U, A, S)
        F = D.normalized_objective(Z.to(t2.dtype), Uk, t2, lam, zref)
        gA, gS = torch.autograd.grad(F.sum(), [A, S])
    gA = 0.5 * (gA - gA.transpose(-1, -2))
    gS = 0.5 * (gS + gS.transpose(-1, -2))
    gn = torch.sqrt((gA.pow(2) + gS.pow(2)).flatten(1).mean(1)).clamp_min(1e-12)[:, None, None, None]
    gA, gS = gA / gn, gS / gn
    occw = U[..., :nocc, :].abs().pow(2).sum(-2)                       # (B,R,n)
    zd = torch.diagonal(Z, dim1=-2, dim2=-1)                            # (B,R,n)
    feats = []
    for k in range(R):
        feats += [G[:, k].real, G[:, k].imag, GR[:, k].real, GR[:, k].imag, gA[:, k].double(), gS[:, k].double(),
                  Z[:, k].double() * 10.0, M[:, k, k].abs(), M[:, k, k].real, M[:, k, k].imag]
    feats += [M[:, 0, 1].real, M[:, 0, 1].imag]
    feats.append(mask.double().expand(B, n, n))
    dist = (torch.arange(n, device=U.device)[:, None] - torch.arange(n, device=U.device)[None, :]).abs().clamp(max=3)
    feats += [(dist == j).double().expand(B, n, n) for j in range(4)]
    pair = torch.stack(feats, -1).float()                               # (B,n,n,N_FEAT_PAIR)
    node = torch.cat([occw.permute(0, 2, 1), zd.permute(0, 2, 1) * 10.0], -1).float()   # (B,n,2R)
    obj = D.normalized_objective(Z, U, t2, lam, zref)
    return pair, node, Z, obj


class SlotNet(nn.Module):
    def __init__(self, d: int = 128, layers: int = 4, heads: int = 8, n_reps: int = 2, kscale: float = 1.0,
                 max_recycles: int = 16, real: bool = False, n_chem_node: int = 0, n_chem_pair: int = 0):
        super().__init__()
        self.R, self.kscale, self.real = n_reps, kscale, real
        node_dim = 2 * n_reps
        self.inp = nn.Sequential(nn.Linear(N_FEAT_PAIR + 2 * node_dim + 2 * N_POS + 2 * n_chem_node + n_chem_pair, d),
                                 nn.GELU(), nn.Linear(d, d))
        self.recycle_norm = nn.LayerNorm(d)
        self.t_emb = nn.Embedding(max_recycles, d)
        cfg = ModelConfig(embed_dim=d, num_layers=layers, num_heads=heads, dropout=0.0, attention_dropout=0.0)
        self.layers = nn.ModuleList([TransformerLayer(cfg) for _ in range(layers)])
        self.norm = nn.LayerNorm(d)

        def head():
            h = nn.Sequential(nn.Linear(d, d // 2), nn.GELU(), nn.Linear(d // 2, 1))
            nn.init.zeros_(h[-1].weight)
            nn.init.zeros_(h[-1].bias)
            return h
        self.f = nn.ModuleList([head() for _ in range(n_reps)])
        self.h = nn.ModuleList([head() for _ in range(n_reps)])

    def step(self, pair, node, prev, t, chem=None):
        B, n, _, _ = pair.shape
        pos = sinusoid(n, N_POS, pair.device)
        xs = [pair, node[:, :, None, :].expand(B, n, n, -1), node[:, None, :, :].expand(B, n, n, -1),
              pos[None, :, None, :].expand(B, n, n, -1), pos[None, None, :, :].expand(B, n, n, -1)]
        if chem is not None:
            cn, cp = chem
            xs += [cn[:, :, None, :].expand(B, n, n, -1), cn[:, None, :, :].expand(B, n, n, -1), cp]
        x = torch.cat(xs, -1)
        z = self.inp(x) + self.t_emb.weight[min(t, self.t_emb.num_embeddings - 1)]
        if prev is not None:
            z = z + self.recycle_norm(prev)
        for L in self.layers:
            z = L(z)
        z = self.norm(z)
        zt = z.transpose(1, 2)
        dK = []
        for k in range(self.R):
            re = torch.tanh(self.f[k](z).squeeze(-1) - self.f[k](zt).squeeze(-1)) * self.kscale
            if self.real:                       # real sector: real antisymmetric generator only
                im = torch.zeros_like(re)
            else:
                im = torch.tanh(self.h[k](z).squeeze(-1) + self.h[k](zt).squeeze(-1)) * self.kscale
            dK.append(torch.complex(re, im))
        return torch.stack(dK, 1), z                                    # (B,R,n,n), state

    def forward(self, U0: Tensor, t2: Tensor, mask: Tensor, lam: float, zref: Tensor, T: int = 1,
                detach_state: bool = False, chem=None, gen_mask: Tensor | None = None):
        """Run T recycles from U0; returns list of (U_t, Z*_t, obj_t) for t = 0..T (t=0 is the start).

        real=True: U0 = Phi B is in the real sector and stays there (real dK, no phase re-canonicalization).
        chem = (node (B,n,Fn), pair (B,n,n,Fp)) static chemistry features; gen_mask (B,n,n) restricts dK."""
        U = U0 if self.real else canon_phases(U0)
        traj = []
        prev = None
        for t in range(T + 1):
            pair, node, Z, obj = features(U, t2, mask, lam, zref)
            traj.append((U, Z, obj))
            if t == T:
                break
            dK, prev = self.step(pair, node, prev, t, chem)
            if gen_mask is not None:
                dK = dK * gen_mask[:, None].to(dK.real.dtype)
            U = U @ torch.linalg.matrix_exp(dK.to(U.dtype))
            if not self.real:
                U = canon_phases(U)
            if detach_state:
                U = U.detach()
        return traj
