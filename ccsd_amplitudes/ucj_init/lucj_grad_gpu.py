"""Exact LUCJ energy AND gradient on one GPU (adjoint method), for direct optimization from a given start.

The state is the project's LUCJ (pretrain.rl.energy.make_ucj_op):
    psi = W_f  D_L W_L ... D_1 W_1 |HF>,   W_1 = U_1^dag, W_r = U_r^dag U_{r-1}, W_f = F U_L,
D_r = exp(i J_r) (diagonal Coulomb, square connectivity: J_aa on (p, p+1), J_ab on (p, p)), F = exp(t1 - t1^dag)
(fixed, not optimized).  Parameters x = ffsim UCJOpSpinBalanced.to_parameters(interaction_pairs=square) without a
final rotation, i.e. per layer: n^2 parameters of log U (ffsim antihermitian layout), n-1 J_aa, n J_ab (task 2's x).

Gradient (adjoint).  E = <psi|H|psi>; lambda_m = H psi; going back through the gates g_j:
  * rotation W (right perturbation W -> W e^{eps A}):  dE = 2 Re sum_pq A_pq rho_pq,
    rho_pq = <lambda_{j-1}| a^dag_p a_q |phi_{j-1}> (spin-summed transition 1-RDM, before the gate);
  * phase exp(i J):  dE/dJ_term = -2 Im <lambda_j| term |phi_j> (term = n_pa n_qa + n_pb n_qb, or n_pa n_pb);
  then phi_{j-1} = g_j^dag phi_j, lambda_{j-1} = g_j^dag lambda_j.  The U-gradient goes through W_1 = U_1^dag,
  W_2 = U_2^dag U_1, W_f = F U_L and U = expm(K) (adjoint Frechet derivative: tr(M L(K, E)) = tr(L(K, M) E)).

H psi on the GPU, HF-shifted like pretrain.rl.gpu_energy.factorize_hamiltonian:
  H = c + E0 + Ha' (x) 1 + 1 (x) Ha' + sum_PQ (P|Q) Ft^a_P Ft^b_Q,
  Ha' = sum_P g_P Ft^a_P + 1/2 sum_PQ (P|Q) Ft^a_P Ft^a_Q   (dense dim x dim, built once on the CPU),
  Ft_P = F_P - <HF|F_P|HF>, packed pairs P = (p <= q), F_pq = E_pq + E_qp.
  With psi = psi^T (spin-balanced closed shell): sigma = (c + E0) psi + Ha' psi + psi Ha'^T + V_ab(psi).
  V_ab: W[b, P, :] = sum_Q (P|Q) (Ft^b_Q psi)[:, b] is the energy engine's X tile with the (P|Q) rows in place of
  the factor rows (one batched GEMM per tile); then V_ab[a, b] = sum_e s(a,e) W[b, P(a,e), src(a,e)] + diagonal.

    eng = LUCJGradGPU.from_npz(name, t1=t1)
    E, g = eng.energy_grad(x)        # x: float64 (n_params,), g: float64
"""
from __future__ import annotations

import math
import time
from functools import lru_cache

import numpy as np
import scipy.linalg
import torch

from pretrain.rl.gpu_energy import (LUCJEnergyGPU, factorize_hamiltonian, givens_clustered, split_segments,
                                    string_tables)


def square_pairs(norb: int):
    return [(p, p + 1) for p in range(norb - 1)], [(p, p) for p in range(norb)]


def n_params(norb: int, n_reps: int) -> int:
    return n_reps * (norb ** 2 + (norb - 1) + norb)


def antiherm_from_params(v: np.ndarray, n: int) -> np.ndarray:
    """ffsim.linalg.util.antihermitians_from_parameters for one matrix."""
    nt = n * (n - 1) // 2
    K = np.zeros((n, n), dtype=complex)
    r, c = np.triu_indices(n)
    K[r, c] = 1j * v[nt:]
    K[c, r] = 1j * v[nt:]
    r, c = np.triu_indices(n, k=1)
    K[r, c] += v[:nt]
    K[c, r] -= v[:nt]
    return K


def grad_params_from_N(N: np.ndarray) -> np.ndarray:
    """dE/dv for dE = 2 Re tr(N dK) and K = antiherm_from_params(v)."""
    n = N.shape[0]
    r, c = np.triu_indices(n, k=1)
    g_re = 2.0 * (N[c, r] - N[r, c]).real
    r2, c2 = np.triu_indices(n)
    s = np.where(r2 == c2, N[r2, c2], N[c2, r2] + N[r2, c2])
    g_im = 2.0 * (1j * s).real
    return np.concatenate([g_re, g_im])


def unpack_x(x: np.ndarray, norb: int, n_reps: int):
    """x -> (K (R,n,n) complex antihermitian, Jaa (R, n-1), Jab (R, n))."""
    n = norb
    Ks, Jaa, Jab = [], [], []
    i = 0
    for _ in range(n_reps):
        Ks.append(antiherm_from_params(x[i:i + n * n], n))
        i += n * n
        Jaa.append(np.asarray(x[i:i + n - 1], dtype=float))
        i += n - 1
        Jab.append(np.asarray(x[i:i + n], dtype=float))
        i += n
    return np.stack(Ks), np.stack(Jaa), np.stack(Jab)


def uz_from_x(x: np.ndarray, norb: int, n_reps: int):
    """x -> (U, Z) in the project convention (Z: J_aa on (p,p+1) and J_ab on the diagonal, one matrix)."""
    K, Jaa, Jab = unpack_x(x, norb, n_reps)
    U = scipy.linalg.expm(K)
    Z = np.zeros((n_reps, norb, norb))
    p = np.arange(norb - 1)
    Z[:, p, p + 1] = Jaa
    Z[:, p + 1, p] = Jaa
    Z[:, np.arange(norb), np.arange(norb)] = Jab
    return U, Z


def x_from_uz(U: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """(U (R,n,n), Z (R,n,n)) -> x (polar-projected U, principal log)."""
    R, n, _ = U.shape
    out = []
    for r in range(R):
        W_, _, Vh = np.linalg.svd(U[r])
        K = scipy.linalg.logm(W_ @ Vh)
        K = 0.5 * (K - K.conj().T)
        nt = n * (n - 1) // 2
        ru, cu = np.triu_indices(n, k=1)
        r2, c2 = np.triu_indices(n)
        out += [K[ru, cu].real, K[r2, c2].imag, Z[r, np.arange(n - 1), np.arange(1, n)], np.diag(Z[r])]
        assert len(out[-4]) == nt
    return np.concatenate(out).astype(np.float64)


@lru_cache(maxsize=None)
def link_tables(norb: int, k: int):
    """Target-indexed single excitations (off-diagonal): for string a and e < k(n-k):
    (E_{p q} phi)[a] += sign * phi[src]  with  p = P1[a,e], q = Q1[a,e]."""
    from pyscf.fci import cistring
    link = cistring.gen_linkstr_index(range(norb), k)
    off = link[:, :, 0] != link[:, :, 1]
    dim = link.shape[0]
    lo = link[off].reshape(dim, -1, 4)
    # pyscf link[a] = (x, y, str1, sign): E_xy|a> = sign|str1>  =>  <a|E_yx|str1> = sign
    P1 = lo[:, :, 1].astype(np.int64)
    Q1 = lo[:, :, 0].astype(np.int64)
    src = lo[:, :, 2].astype(np.int64)
    sgn = lo[:, :, 3].astype(np.float64)
    return P1, Q1, src, sgn


class LUCJGradGPU:
    def __init__(self, one_body, two_body, constant, norb, nelec, t1=None, n_reps=2, device="cuda",
                 e_hf=None, e_ccsd=None, name=None, wtile=None, verbose=False):
        self.norb, self.k = int(norb), int(nelec[0])
        assert nelec[0] == nelec[1]
        self.n_reps = n_reps
        self.device = torch.device(device)
        if self.device.index is None:
            self.device = torch.device("cuda", torch.cuda.current_device())
        self.e_hf, self.e_ccsd, self.name = e_hf, e_ccsd, name
        self.cdtype = torch.complex64
        n, k = self.norb, self.k
        t0 = time.time()
        with torch.cuda.device(self.device):
            # kernels / tables of the energy engine (no state buffer, tiny workspace)
            self.eng = LUCJEnergyGPU(one_body, two_body, constant, n, nelec, device=self.device,
                                     alloc_state=False, tile=32)
            self.be = self.eng.backend
            assert self.be.name == "numba", "the numba backend is required"
            self.g = self.eng.g
            self.dim = dim = self.eng.dim
        if t1 is not None:
            from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
            self.F = orbital_rotation_from_t1_amplitudes(np.asarray(t1, dtype=np.float64))
        else:
            self.F = np.eye(n)
        self._build_hamiltonian(np.asarray(one_body, np.float64), np.asarray(two_body, np.float64),
                                float(constant))
        self.n_params = n_params(n, n_reps)
        self.n_calls = 0
        self.t_build = time.time() - t0
        self.verbose = verbose
        dev = self.device
        with torch.cuda.device(dev):
            self.phi = torch.empty((dim, dim), device=dev, dtype=self.cdtype)
            self.lam = torch.empty((dim, dim), device=dev, dtype=self.cdtype)
            self.sig = torch.empty((dim, dim), device=dev, dtype=self.cdtype)
            free = torch.cuda.mem_get_info(dev)[0]
            K = self.K2
            NP = self.npair
            # W tile (B, NP+1, dim) + gather G (B, K, dim) + gathered (B, dim, K-1)  complex
            per_b = 8 * dim * (NP + 1 + K + K)
            self.B = wtile or max(8, min(256, int(0.5 * (free - 0.5e9) / per_b)))
        if verbose:
            print(f"[grad] {name} norb {n} k {k} dim {dim} npair {self.npair} B {self.B} build {self.t_build:.1f}s",
                  flush=True)

    @classmethod
    def from_npz(cls, name, ham_dir=None, **kw):
        from pathlib import Path
        p = Path(ham_dir or (Path(__file__).resolve().parents[1] / "rhf_hamiltonians")) / f"{name}.npz"
        d = np.load(p)
        return cls(d["one_body"], d["two_body"], float(d["constant"]), int(d["norb"]),
                   (int(d["nelec_a"]), int(d["nelec_b"])), e_hf=float(d["e_hf"]), e_ccsd=float(d["e_ccsd"]),
                   name=Path(p).stem, **kw)

    # ------------------------------------------------------------------ Hamiltonian
    def _build_hamiltonian(self, h1, eri, const):
        n, k, dev = self.norb, self.k, self.device
        tb = string_tables(n, k)
        dim, K = tb["dim"], tb["K"]
        self.K2 = K
        _, _, _, E0, _ = factorize_hamiltonian(h1, eri, k)
        iu = np.triu_indices(n)
        npair = len(iu[0])
        self.npair = npair
        PQ = eri[iu[0], iu[1]][:, iu[0], iu[1]]
        PQ = 0.5 * (PQ + PQ.T)
        hp = h1 - 0.5 * np.einsum("prrq->pq", eri)
        nP = np.zeros(npair)
        nP[(iu[0] == iu[1]) & (iu[0] < k)] = 2.0
        gP = hp[iu] + PQ @ nP
        self.E0, self.const = E0, const
        hf = np.zeros(n)
        hf[:k] = 1.0
        occ = tb["occ"].astype(np.float64)
        src, exc_p, exc_s, diagP = tb["src"], tb["exc_p"], tb["exc_s"], tb["diagP"]   # src[:, -1] = self
        # C2[a, P, e]: sum_Q (P|Q) (Ft_Q phi)[a] = sum_e C2[a,P,e] phi[src(a,e)];  extra row P = npair: g
        Vext = np.concatenate([PQ, gP[None]], 0)                                # (npair+1, npair)
        C2 = np.empty((dim, npair + 1, K), dtype=np.float64)
        C2[:, :, :-1] = np.transpose(Vext[:, exc_p], (1, 0, 2)) * exc_s[:, None, :]
        C2[:, :, -1] = (occ - hf) @ Vext[:, diagP].T
        # dense Ha' = Bt_g + 1/2 sum_P Ft_P Y_P,  Y_P[s, src(s,e')] = C2[s, P, e']
        #   (Ft_P Y_P)[a, :] = sum_e s(a,e) Y_{P(a,e)}[src(a,e), :]  +  sum_q (occ[a,q]-hf[q]) Y_{qq}[a, :]
        Ha = np.zeros(dim * dim)
        rows = np.repeat(np.arange(dim), K)
        Ha += np.bincount(rows * dim + src.reshape(-1), weights=C2[:, npair, :].reshape(-1), minlength=dim * dim)
        Koff = K - 1
        chunk = max(1, int(4e7 // (Koff * K)))
        for a0 in range(0, dim, chunk):
            a1 = min(dim, a0 + chunk)
            s_ae = src[a0:a1, :Koff]                                           # (A, Koff)
            P_ae = exc_p[a0:a1]                                                # (A, Koff)
            w = exc_s[a0:a1][:, :, None] * C2[s_ae, P_ae, :]                   # (A, Koff, K)
            cols = src[s_ae]                                                   # (A, Koff, K)
            rr = np.arange(a0, a1)[:, None, None]
            Ha += 0.5 * np.bincount((rr * dim + cols).reshape(-1), weights=w.reshape(-1), minlength=dim * dim)
            # diagonal part: sum_q (occ[a,q]-hf[q]) Y_qq[a, :]
            dq = (occ[a0:a1] - hf)                                             # (A, n)
            wd = np.einsum("aq,aqe->ae", dq, C2[a0:a1][:, diagP, :])           # (A, K)
            Ha += 0.5 * np.bincount((np.arange(a0, a1)[:, None] * dim + src[a0:a1]).reshape(-1),
                                    weights=wd.reshape(-1), minlength=dim * dim)
        Ha = Ha.reshape(dim, dim)
        self.Ha_asym = float(np.abs(Ha - Ha.T).max())
        Ha = 0.5 * (Ha + Ha.T)
        with torch.cuda.device(dev):
            self.Ha = torch.as_tensor(Ha, device=dev, dtype=torch.float32)
            self.C2 = torch.as_tensor(C2[:, :npair, :], device=dev, dtype=torch.float32).contiguous()
            # gather index into W (B, npair*dim) for the alpha application: P(a,e)*dim + src(a,e)
            self.widx = (torch.as_tensor(exc_p, device=dev) * dim
                         + torch.as_tensor(src[:, :Koff], device=dev)).reshape(-1)
            self.wsgn = torch.as_tensor(exc_s, device=dev, dtype=torch.float32)
            # the same alpha application as one sparse matrix S (dim, npair*dim):
            #   S[a, P(a,e)*dim + src(a,e)] = s(a,e);  S[a, (qq)*dim + a] = occ[a,q] - hf[q]
            rr = np.repeat(np.arange(dim), Koff)
            cc = (exc_p * dim + src[:, :Koff]).reshape(-1)
            vv = exc_s.reshape(-1)
            aq = np.nonzero(occ - hf)
            rr = np.concatenate([rr, aq[0]])
            cc = np.concatenate([cc, diagP[aq[1]] * dim + aq[0]])
            vv = np.concatenate([vv, (occ - hf)[aq]])
            Ssp = torch.sparse_coo_tensor(np.stack([rr, cc]), vv, (dim, npair * dim)).coalesce()
            self.S = Ssp.to(device=dev, dtype=torch.float32).to_sparse_csr()
            self.dq = torch.as_tensor(occ - hf, device=dev, dtype=torch.float32)   # (dim, n)
            self.diagP_t = torch.as_tensor(diagP, device=dev)
            self.src_t = torch.as_tensor(src, device=dev)
            P1, Q1, lsrc, lsgn = link_tables(n, k)
            self.lP = torch.as_tensor(P1, device=dev)
            self.lQ = torch.as_tensor(Q1, device=dev)
            self.lsrc = torch.as_tensor(lsrc, device=dev)
            self.lsgn = torch.as_tensor(lsgn, device=dev, dtype=torch.float32)
            self.occf = torch.as_tensor(occ, device=dev, dtype=torch.float32)
            self.occd = torch.as_tensor(occ, device=dev, dtype=torch.float64)
        del Ha, C2

    def sigma(self, psi, out):
        """out <- (H - c - E0) psi for a symmetric CI matrix psi (complex64)."""
        dim, B, npair = self.dim, self.B, self.npair
        Kfull = self.K2
        prev = torch.backends.cuda.matmul.allow_tf32
        torch.backends.cuda.matmul.allow_tf32 = False
        try:
            pr = torch.view_as_real(psi).view(dim, 2 * dim)
            orr = torch.view_as_real(out).view(dim, 2 * dim)
            torch.mm(self.Ha, pr, out=orr)                                     # Ha' psi
            # + psi Ha'^T = (Ha' psi)^T  since psi = psi^T and Ha' = Ha'^T
            out.add_(out.T.clone())
            for b0 in range(0, dim, B):
                b1 = min(dim, b0 + B)
                nb = b1 - b0
                idx = self.src_t[b0:b1].reshape(-1)
                G = torch.index_select(psi, 0, idx)                            # (nb*K, dim)
                Gr = torch.view_as_real(G).view(nb, Kfull, 2 * dim)
                Wr = torch.bmm(self.C2[b0:b1], Gr)                             # (nb, npair, 2dim)
                Wt = Wr.view(nb, npair, dim, 2).permute(1, 2, 0, 3).reshape(npair * dim, 2 * nb)
                v = torch.sparse.mm(self.S, Wt)                                # (dim, 2nb): [a, (b, re/im)]
                out[:, b0:b1] += torch.view_as_complex(v.view(dim, nb, 2))
                del G, Gr, Wr, Wt, v
        finally:
            torch.backends.cuda.matmul.allow_tf32 = prev
        return out

    # ------------------------------------------------------------------ gates on a buffer
    def _rotate(self, X, W):
        """X <- Lambda(W) X Lambda(W)^T  (symmetric X)."""
        app, D = givens_clustered(np.asarray(W, dtype=np.complex128), self.eng.t_top)
        segs = split_segments(app, self.norb, self.eng.t_top)
        dA = torch.exp(self.occd.to(torch.complex128)
                       @ torch.as_tensor(np.log(D.astype(np.complex128)), device=self.device))
        ones = torch.ones(self.norb, device=self.device, dtype=torch.complex128)
        Xb = self.be.arr(X)
        self.be.phase_diag(Xb, self.g["strs_nb"], dA.contiguous(), ones, self.norb, False)
        if self.eng.chunks is not None:
            useg = self.be.upload_segments(segs)
            self.be.givens_minor(Xb, self.dim, useg, self.eng.pairs_b, self.eng.chunks)
            self.be.transpose(Xb, self.dim)
            self.be.givens_minor(Xb, self.dim, useg, self.eng.pairs_b, self.eng.chunks)
        else:
            rots = [r for _, seg in segs for r in seg]
            self.be.givens_sequence(Xb, self.dim, rots, self.eng.pairs_b)
            self.be.transpose(Xb, self.dim)
            self.be.givens_sequence(Xb, self.dim, rots, self.eng.pairs_b)

    def _phase(self, X, jaa, jab, sign=1.0):
        """X[a,b] *= exp(i sign [phi(a) + phi(b) + sum_p jab_p n_p(a) n_p(b)]),  phi(s) = sum_p jaa_p n_p n_{p+1}."""
        occ = self.occd
        phi = (occ[:, :-1] * occ[:, 1:]) @ torch.as_tensor(jaa, device=self.device, dtype=torch.float64)
        rowfac = torch.polar(torch.ones_like(phi), sign * phi)
        ez = torch.polar(torch.ones(self.norb, device=self.device, dtype=torch.float64),
                         sign * torch.as_tensor(jab, device=self.device, dtype=torch.float64))
        self.be.phase_diag(self.be.arr(X), self.g["strs_nb"], rowfac.contiguous(), ez.contiguous(), self.norb,
                           False)

    def _slater(self, X, W1):
        occl = self.g["occ_list"]
        W1t = torch.as_tensor(np.asarray(W1, dtype=np.complex128), device=self.device)
        v = torch.linalg.det(W1t[occl][:, :, : self.k])
        X.copy_((v[:, None] * v[None, :]).to(self.cdtype))

    # ------------------------------------------------------------------ RDM-type contractions
    def _trdm1(self, lam, phi):
        """rho_pq = <lam| sum_s a^dag_ps a_qs |phi> for symmetric lam, phi (= 2x the alpha part)."""
        n, dim = self.norb, self.dim
        rho = torch.zeros(n * n, device=self.device, dtype=torch.complex128)
        T = max(1, int(2e8 // (8 * dim * self.lsrc.shape[1])))
        lc = lam.conj()
        for a0 in range(0, dim, T):
            a1 = min(dim, a0 + T)
            Gs = phi[self.lsrc[a0:a1]]                                         # (A, Koff, dim)
            S = torch.einsum("ab,aeb->ae", lc[a0:a1], Gs) * self.lsgn[a0:a1]    # (A, Koff)
            ii = (self.lP[a0:a1] * n + self.lQ[a0:a1]).reshape(-1)
            rho.index_add_(0, ii, S.reshape(-1).to(torch.complex128))
            del Gs, S
        r = (lc * phi).sum(1).to(torch.complex128)                            # (dim,)
        rho = rho.view(n, n)
        rho += torch.diag(self.occd.T.to(torch.complex128) @ r)
        return 2.0 * rho.cpu().numpy()

    def _phase_grads(self, lam, phi):
        """(-2 Im <lam| n_pa n_{p+1,a} + n_pb n_{p+1,b} |phi>,  -2 Im <lam| n_pa n_pb |phi>)."""
        w = (lam.conj() * phi)                                                # (dim, dim)
        r = w.sum(1).to(torch.complex128)
        occ = self.occd.to(torch.complex128)
        same = 2.0 * ((occ[:, :-1] * occ[:, 1:]).T @ r)                        # alpha + beta
        oc = self.occf.to(w.dtype)
        ab = ((oc.T @ w) * oc.T).sum(1).to(torch.complex128)                   # sum_ab w[a,b] n_p(a) n_p(b)
        return (-2.0 * same.imag).cpu().numpy(), (-2.0 * ab.imag).cpu().numpy()

    # ------------------------------------------------------------------ public
    @torch.no_grad()
    def energy_grad(self, x, want_grad=True):
        n, R = self.norb, self.n_reps
        x = np.asarray(x, dtype=np.float64)
        Ks, Jaa, Jab = unpack_x(x, n, R)
        Us = [scipy.linalg.expm(K) for K in Ks]
        Ws = [Us[0].conj().T] + [Us[r].conj().T @ Us[r - 1] for r in range(1, R)] + [self.F @ Us[-1]]
        phi, lam, sig = self.phi, self.lam, self.sig
        with torch.cuda.device(self.device):
            self._slater(phi, Ws[0])
            for r in range(R):
                self._phase(phi, Jaa[r], Jab[r])
                self._rotate(phi, Ws[r + 1])
            self.be.symmetrize(self.be.arr(phi), self.dim)
            self.sigma(phi, sig)
            nrm = float(torch.vdot(phi.reshape(-1), phi.reshape(-1)).real)
            e_rel = float(torch.vdot(phi.reshape(-1), sig.reshape(-1)).real) / nrm
            E = self.const + self.E0 + e_rel
            self.n_calls += 1
            if not want_grad:
                return E
            # adjoint sweep; lambda = (H - E) psi (the shift removes nothing from the exact gradient but
            # keeps it orthogonal to psi under float32 rounding)
            lam.copy_(sig)
            lam.add_(phi, alpha=-e_rel)
            grads_W = [None] * (R + 1)
            gJaa = np.zeros((R, n - 1))
            gJab = np.zeros((R, n))
            for r in range(R, 0, -1):
                # rotation W_{r} (index r in Ws; r = R is W_f)
                Wd = Ws[r].conj().T
                self._rotate(phi, Wd)
                self._rotate(lam, Wd)
                grads_W[r] = self._trdm1(lam, phi)
                # phase layer r-1 (0-based): after-gate quantities
                gaa, gab = self._phase_grads(lam, phi)
                gJaa[r - 1] = gaa
                gJab[r - 1] = gab
                self._phase(phi, Jaa[r - 1], Jab[r - 1], sign=-1.0)
                self._phase(lam, Jaa[r - 1], Jab[r - 1], sign=-1.0)
            # first rotation acts on HF
            self._rotate(lam, Ws[0].conj().T)
            hf = torch.zeros_like(phi)
            hf[0, 0] = 1.0
            grads_W[0] = self._trdm1(lam, hf)
        # chain rule  dE = 2 Re sum_j tr(dW_j Gam_j),  Gam_j = rho_j^T W_j^dag
        Gam = [grads_W[j].T @ Ws[j].conj().T for j in range(R + 1)]
        M = [np.zeros((n, n), dtype=complex) for _ in range(R)]
        # W_0 = U_0^dag:      Re tr(dU^dag G) = Re tr(G^dag dU)
        M[0] += Gam[0].conj().T
        for r in range(1, R):
            # W_r = U_r^dag U_{r-1}:  tr(dU_r^dag U_{r-1} G) + tr(U_r^dag dU_{r-1} G)
            M[r] += (Us[r - 1] @ Gam[r]).conj().T
            M[r - 1] += Gam[r] @ Us[r].conj().T
        M[R - 1] += Gam[R] @ self.F                                           # W_f = F U_L
        grad = []
        for r in range(R):
            N = scipy.linalg.expm_frechet(Ks[r], M[r], compute_expm=False)
            grad += [grad_params_from_N(N), gJaa[r], gJab[r]]
        return E, np.concatenate(grad)

    def corr_pct(self, E):
        return 100.0 * (self.e_hf - E) / (self.e_hf - self.e_ccsd)
