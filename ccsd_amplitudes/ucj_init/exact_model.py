"""Exact-energy drop-in for Avijit's `_LocalModel` (small systems, dense state vector, torch autograd).

The proposed fix keeps RUCJ.kernel (trust-region Newton, Steihaug-CG, accept/reject, convergence tests) and
the parameterization (accumulated orbital rotation U, occupied-virtual kappa steps, packed symmetric g) unchanged,
and replaces only the BCH2 energy model by the exact energy of the same state,

    E(kappa, dg) = <HF| e^{-iX} H e^{iX} |HF>,   X = R (G(g + dg)),   R = U expm(K(kappa)),

with exact gradient and Hessian-vector products at kappa = 0 (autograd, complex128).  The state is built as in
ffsim (verified in bch2_reproduce.py): psi = Lam(R^T) exp(i J) Lam(R^T)^dag |HF>, J_aa = J_ab = g, where
Lam(W) is the action of the orbital rotation W on k-electron strings: Lam(W)[a, b] = det W[occ(a), occ(b)],
and Lam(R^T) = expm(-dLam(K)) Lam(U^T) with dLam(K) = sum_pq K_pq E_pq.

    mod = load_avijit(); install(mod, mf)   # then mod.UCJ(mf).kernel(g0, u0) optimizes the exact energy
"""
from __future__ import annotations

import itertools

import numpy as np
import torch

torch.set_default_dtype(torch.float64)


class DenseFCI:
    """String-space tables and the Hamiltonian action for a small closed-shell system (complex128, CPU)."""

    def __init__(self, h1, eri, const, n, k):
        from pyscf.fci import cistring
        self.n, self.k = n, k
        strs = cistring.make_strings(range(n), k)
        self.dim = len(strs)
        occ = np.array([[(s >> p) & 1 for p in range(n)] for s in strs], dtype=float)
        self.occ_np = occ
        self.occ_idx = [np.nonzero(o)[0] for o in occ]
        self.occ = torch.as_tensor(occ)
        # E_pq on alpha strings: E[p, q][a, b] = <a| a^dag_p a_q |b>
        E = np.zeros((n, n, self.dim, self.dim))
        link = cistring.gen_linkstr_index(range(n), k)
        for b in range(self.dim):
            for (p, q, a, sgn) in link[b]:
                E[p, q, a, b] += sgn
        self.E = torch.as_tensor(E)
        hp = h1 - 0.5 * np.einsum("prrq->pq", eri)
        Ha = np.einsum("pq,pqab->ab", hp, E) + 0.5 * np.einsum("pqrs,pqac,rscb->ab", eri, E, E, optimize=True)
        self.Ha = torch.as_tensor(Ha).to(torch.complex128)
        self.eri = torch.as_tensor(eri.reshape(n * n, n * n)).to(torch.complex128)
        self.Ec = self.E.reshape(n * n, self.dim, self.dim).to(torch.complex128)
        self.const = const
        hf = torch.zeros(self.dim, self.dim, dtype=torch.complex128)
        hf[0, 0] = 1.0
        self.hf = hf

    def sigma(self, psi):
        Y = torch.einsum("ab,rcb->rac", psi, self.Ec)                  # psi E_rs^T
        Z = torch.einsum("pr,rac->pac", self.eri, Y)
        V = torch.einsum("pab,pbc->ac", self.Ec, Z)
        return self.Ha @ psi + psi @ self.Ha.T + V + self.const * psi

    def energy(self, psi):
        return torch.vdot(psi.reshape(-1), self.sigma(psi).reshape(-1)).real / torch.vdot(
            psi.reshape(-1), psi.reshape(-1)).real

    def lam(self, W):
        """Lam(W) by determinants (no autograd), W real or complex (n, n)."""
        W = np.asarray(W)
        idx = self.occ_idx
        out = np.empty((self.dim, self.dim), dtype=complex)
        for a in range(self.dim):
            sub = W[idx[a]]
            out[a] = np.linalg.det(np.stack([sub[:, idx[b]] for b in range(self.dim)]))
        return torch.as_tensor(out)

    def gen(self, K):
        """dLam(K) = sum_pq K_pq E_pq (torch, differentiable)."""
        return torch.einsum("pq,pqab->ab", K.to(torch.complex128), self.E.to(torch.complex128))

    def phase(self, g):
        occ = self.occ
        same = 0.5 * ((occ @ g) * occ).sum(1)            # 1/2 sum_pq g n_p n_q  (incl. 1/2 g_pp n_p)
        return same[:, None] + same[None, :] + occ @ g @ occ.T


def ucj_state(fci: DenseFCI, A0, K, g):
    """psi = Lam(R^T) exp(i J(g)) Lam(R^T)^dag HF,  Lam(R^T) = expm(-dLam(K)) A0 (A0 = Lam(U^T))."""
    A = torch.linalg.matrix_exp(-fci.gen(K)) @ A0
    phi = A.conj().T @ fci.hf @ A.conj()
    phi = phi * torch.exp(1j * fci.phase(g))
    return A @ phi @ A.T


def install(mod, fci: DenseFCI):
    """Replace mod._LocalModel by the exact-energy model (same interface as used by RUCJ.kernel)."""
    base = mod._LocalModel

    class ExactLocalModel:
        def __init__(self, parent, g, rotation):
            self.parent, self.g, self.rotation = parent, np.asarray(g), np.asarray(rotation)
            n = parent.nmo
            self.A0 = fci.lam(self.rotation.T)
            self.g_t = torch.as_tensor(self.g)
            self.orbital_mask = np.ones(parent.nkappa, dtype=bool)
            iu = parent._g_indices
            self.iu = (torch.as_tensor(iu[0]), torch.as_tensor(iu[1]))
            self.n = n
            x0 = torch.zeros(parent.nparam, requires_grad=True)
            E = self._energy(x0)
            (gr,) = torch.autograd.grad(E, x0, create_graph=True)
            self._x0, self._E, self._gr = x0, float(E), gr
            self._grad = self.project(gr.detach().numpy())

        def _energy(self, x):
            p = self.parent
            pk, pg = x[:p.nkappa], x[p.nkappa:]
            K = torch.zeros(self.n, self.n)
            blk = pk.reshape(p.nvir, p.nocc)
            K = K.index_put((torch.arange(p.nocc, self.n)[:, None], torch.arange(p.nocc)[None, :]), blk)
            K = K - K.T
            dg = torch.zeros(self.n, self.n).index_put(self.iu, pg)
            dg = dg + dg.T - torch.diag(torch.diagonal(dg))
            psi = ucj_state(fci, self.A0, K, self.g_t + dg)
            return fci.energy(psi)

        @property
        def energy_elec(self):
            return self._E - self.parent.e_nuc

        project = base.project

        def gradient(self):
            return self._grad.copy()

        def hessp(self, direction):
            d = torch.as_tensor(self.project(direction))
            (hv,) = torch.autograd.grad(self._gr @ d, self._x0, retain_graph=True)
            return self.project(hv.numpy())

    mod._LocalModel = ExactLocalModel
    return ExactLocalModel


def fci_for_mf(mf):
    from ucj_init.common import mo_integrals
    h1, eri = mo_integrals(mf)
    n = h1.shape[0]
    return DenseFCI(h1, eri, mf.mol.energy_nuc(), n, mf.mol.nelectron // 2)
