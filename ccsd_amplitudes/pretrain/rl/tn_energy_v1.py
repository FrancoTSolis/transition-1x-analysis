# Frozen copy of the first tn_energy.py (zip-up MPO application; md5 a2ec0032 before this line). Every n29 TN number
# in docs/rl_larger_molecules.md (Oct 3 2026) was computed with it; tn_energy.py was later replaced by the
# density-matrix engine. Select it in the reward queue with `worker --kind tn --tn-impl v1`.
"""Scalable LUCJ energies with a particle-number-symmetric MPS -- reward engine for RL beyond norb ~17.

Public API
    ev = LUCJEnergyTN(one_body, two_body, constant, norb, nelec, max_bond=64, device="cuda", name=name)
    # basis: "boys" when name is given (geometry), else "er" (integrals only); or an explicit (S, occ_mask)
    E, info = ev.energy(U, Z, t1)
returns <psi|H|psi> (Hartree, incl. constant) of the ffsim state of pretrain/rl/energy.py
    make_ucj_op(Z, U, "square", t1)  (U re-unitarized by its polar factor, as in the energy jobs)
and info (discarded weight sum/max over the optimal truncations, max bond, timings).

Method (what works)
  * Orbital basis of the MPS: occupied and virtual MOs localized SEPARATELY (Edmiston-Ruedenberg from the integrals,
    or Boys/PM from the geometry), sites ordered by a Fiedler vector.  |HF> is then an exact product state and the
    LUCJ correlation is local: at norb 15 the final state needs chi(1e-5 discarded per cut) = 159, against 1123 in
    the MO energy order and >= 393 already for |HF> in the network's chemistry-frame chain order.
  * In that basis S the state is  phi_S = exp(i J_1(n^{V_1})) exp(i J_0(n^{V_0})) |HF_S>,  V_k = S^T U_k: every
    square-mask term is a commuting factor exp(i z n_A n_B) = 1 + (e^{iz}-1) b+_A b+_B b_B b_A of two rotated modes,
    applied EXACTLY as a bond-17 MPO (Jordan-Wigner strings, U(1)xU(1) labels) by zip-up (bond <= zip_margin*chi,
    no relative cutoff: the zip-up Gram is not a Schmidt spectrum) followed by an optimal canonical compression to
    chi.  No orbital rotation of the MPS is ever performed.
  * The t1 rotation and S are absorbed into the Hamiltonian (integrals rotated by final(t1) @ S, real); <H> is the
    block2 quantum-chemistry MPO expectation (SZ, complex MPS, converter SymMPS -> block2 verified to 1e-13).
  * Engines: NpSymMPS (numpy complex128, CPU, low overhead) and GpuSymMPS (torch complex64 tensors on CUDA, labels
    on the host, small sector eigh batched on the CPU).  Both are exact (1e-12) without truncation on random tests.

What does not work (kept for the record): TEBD of the Givens-decomposed circuit in the frame chain order
(SymMPS.apply_orbital_rotation / LUCJEnergyTNGivens).  The Clements networks pass through volume-law intermediate
states (norb 15: 0.29 discarded weight for one rotation at chi 128), and the frame-0 -> frame-1 rotation W has
O(1) elements between orbitals 14+ sites apart in the chain at norb 29, so no small-angle/banded path exists.
Fishman-White preparation (fishman_white_sequence) fixes the determinant part but not W.

See pretrain/rl/tests/ (test_tn_small.py, tn_bench_split.py, tn_rank_eval.py) and results/ for the numbers.
"""
from __future__ import annotations

import math
import time

import numpy as np
import torch

QB = 128                                    # label code = N_alpha * QB + N_beta
SITE_Q = (0, QB, 1, QB + 1)                 # |0>, |a>, |b>, |ab>


# ----------------------------------------------------------------------------------------------- LUCJ structure

def polar_unitary(U: np.ndarray) -> np.ndarray:
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def lucj_layers(U: np.ndarray, Z: np.ndarray, t1: np.ndarray | None = None, *, unitarize: bool = True):
    """(U, Z, t1) -> (rots, dcs, F): orbital rotations and square diagonal-Coulomb layers in application order,
    and the final rotation F (absorbed into the Hamiltonian).  Mirrors ffsim UCJOpSpinBalanced._apply_unitary_."""
    U = np.asarray(U, dtype=np.complex128)
    Z = np.asarray(Z, dtype=np.float64)
    if unitarize:
        U = np.stack([polar_unitary(u) for u in U])
    n_reps, n, _ = U.shape
    cur = np.eye(n, dtype=np.complex128)
    rots, dcs = [], []
    for k in range(n_reps):
        rots.append(U[k].conj().T @ cur)
        Zk = 0.5 * (Z[k] + Z[k].T)
        dcs.append((np.diagonal(Zk, 1).copy(), np.diagonal(Zk).copy()))   # (z_nn (n-1,), z_onsite (n,))
        cur = U[k]
    if t1 is not None:
        from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
        F = orbital_rotation_from_t1_amplitudes(np.asarray(t1, dtype=np.float64)) @ cur
    else:
        F = cur
    return rots, dcs, F


def rotate_hamiltonian(one_body: np.ndarray, two_body: np.ndarray, F: np.ndarray):
    """Integrals of O(F)^dag H O(F): h' = F^dag h F, (pq|rs)' = sum conj(F_ap) F_bq conj(F_cr) F_ds (ab|cd)."""
    F = np.asarray(F, dtype=np.complex128)
    h = F.conj().T @ one_body @ F
    g = np.einsum("abcd,ap->pbcd", two_body, F.conj(), optimize=True)
    g = np.einsum("pbcd,bq->pqcd", g, F, optimize=True)
    g = np.einsum("pqcd,cr->pqrd", g, F.conj(), optimize=True)
    g = np.einsum("pqrd,ds->pqrs", g, F, optimize=True)
    return h, g


# ------------------------------------------------------------------------------------------- local gate matrices

def _creation_ops_4():
    """c^dag_k (16x16) for the 4 modes (i a, i b, i+1 a, i+1 b); basis index = s1*4 + s2, s = n_a + 2 n_b."""
    ops = []
    for k in range(4):
        M = np.zeros((16, 16))
        for idx in range(16):
            s1, s2 = divmod(idx, 4)
            occ = [s1 & 1, (s1 >> 1) & 1, s2 & 1, (s2 >> 1) & 1]
            if occ[k]:
                continue
            sign = (-1) ** sum(occ[:k])
            occ[k] = 1
            new = (occ[0] + 2 * occ[1]) * 4 + (occ[2] + 2 * occ[3])
            M[new, idx] = sign
        ops.append(M)
    return ops


_CDAG4 = _creation_ops_4()


def givens_gate16(g: np.ndarray) -> np.ndarray:
    """16x16 Fock-space matrix (both spins) of the orbital rotation acting as g (2x2) on orbitals (lo, hi):
    c^dag_lo -> g00 c^dag_lo + g10 c^dag_hi,  c^dag_hi -> g01 c^dag_lo + g11 c^dag_hi."""
    c = _CDAG4
    cp = [g[0, 0] * c[0] + g[1, 0] * c[2], g[0, 0] * c[1] + g[1, 0] * c[3],
          g[0, 1] * c[0] + g[1, 1] * c[2], g[0, 1] * c[1] + g[1, 1] * c[3]]
    G = np.zeros((16, 16), dtype=np.complex128)
    vac = np.zeros(16, dtype=np.complex128)
    vac[0] = 1.0
    for idx in range(16):
        s1, s2 = divmod(idx, 4)
        occ = [s1 & 1, (s1 >> 1) & 1, s2 & 1, (s2 >> 1) & 1]
        v = vac
        for k in (3, 2, 1, 0):
            if occ[k]:
                v = cp[k] @ v
        G[:, idx] = v
    return G


def nn_dc_diag16(z: float) -> np.ndarray:
    """diag of exp(i z (n_ia n_ja + n_ib n_jb)) on two adjacent sites."""
    d = np.empty(16, dtype=np.complex128)
    for idx in range(16):
        s1, s2 = divmod(idx, 4)
        d[idx] = np.exp(1j * z * ((s1 & 1) * (s2 & 1) + ((s1 >> 1) & 1) * ((s2 >> 1) & 1)))
    return d


def givens_sequence(R: np.ndarray, tol: float = 1e-12):
    """Orbital rotation R -> ([(lo, g 2x2)] in application order, phases (n,)), O(R) = O(D) O(M_L) ... O(M_1)."""
    from ffsim.linalg import givens_decomposition
    rotations, phases = givens_decomposition(np.asarray(R, dtype=np.complex128), tol=tol)
    seq = []
    for c, s, i, j in rotations:
        lo, hi = (i, j) if i < j else (j, i)
        M = np.eye(2, dtype=np.complex128)
        # G* restricted to (i, j): [[c, conj(s)], [-s, c]] in the (i, j) ordering
        Mij = np.array([[c, np.conj(s)], [-s, c]], dtype=np.complex128)
        if i < j:
            M = Mij
        else:
            M = Mij[::-1, ::-1]
        seq.append((lo, M))
    return seq, np.asarray(phases, dtype=np.complex128)


def slater_sequence(R: np.ndarray, nocc: int, tol: float = 1e-12):
    """Givens sequence preparing O(R)|HF> (occupied orbitals = columns :nocc of R) from |1..1 0..0>, up to a phase."""
    from ffsim.linalg import givens_decomposition_slater
    rotations = givens_decomposition_slater(np.asarray(R[:, :nocc].T, dtype=np.complex128), tol=tol)
    seq = []
    for c, s, i, j in rotations:
        lo, hi = (i, j) if i < j else (j, i)
        Mij = np.array([[c, np.conj(s)], [-s, c]], dtype=np.complex128)
        seq.append((lo, Mij if i < j else Mij[::-1, ::-1]))
    return seq


def fishman_white_sequence(C: np.ndarray, window: int = 12, tol: float = 1e-10):
    """Locality-preserving preparation of the Slater determinant with occupied orbitals = columns of C (n x m).

    Fishman & White, PRB 92, 075132 (2015): sweep i = 0..n-2, diagonalize the correlation matrix on the window
    [i, i+w), rotate its purest eigenvector (eigenvalue closest to 0 or 1) onto site i with a Givens staircase.
    Returns (occ (n,), seq) where seq = [(lo, g)] prepares the determinant from the product state |occ> (up to a
    global phase), and the infidelity 1 - |<D|D_fw>|^2 computed exactly from the orbital overlaps."""
    C = np.asarray(C, dtype=np.complex128)
    n, m = C.shape
    G = C @ C.conj().T                            # projector on the occupied space; O(V): G -> V G V^dag
    recs = []
    occ = np.zeros(n, dtype=int)
    for i in range(n - 1):
        hi = min(n, i + window)
        sub = G[i:hi, i:hi]
        lam, vec = np.linalg.eigh(0.5 * (sub + sub.conj().T))
        purity = np.minimum(lam, 1.0 - lam)
        j = int(np.argmin(purity))
        v = vec[:, j].copy()
        occ[i] = int(lam[j] > 0.5)
        for k in range(hi - i - 1, 0, -1):        # zero v[k] using v[k-1]: rotation on sites (i+k-1, i+k)
            a, b = v[k - 1], v[k]
            r = math.hypot(abs(a), abs(b))
            if r < 1e-300 or abs(b) < tol * max(r, 1e-300):
                continue
            # g (2x2 on (lo, lo+1)) with g @ [a, b] = [r', 0]: rows of a unitary
            g = np.array([[np.conj(a), np.conj(b)], [-b, a]], dtype=np.complex128) / r
            v[k - 1], v[k] = r, 0.0
            lo = i + k - 1
            G[[lo, lo + 1], :] = g @ G[[lo, lo + 1], :]
            G[:, [lo, lo + 1]] = G[:, [lo, lo + 1]] @ g.conj().T
            recs.append((lo, g))
    occ[n - 1] = int(G[n - 1, n - 1].real > 0.5)
    # V = g_K ... g_1 (as n x n), V G0 V^dag ~ diag(occ);  |D> ~ O(V^dag)|occ>: apply g_K^dag first ... g_1^dag last
    seq = [(lo, g.conj().T) for lo, g in reversed(recs)]
    # exact fidelity: prepared orbitals = V^dag e_occ
    Vd = np.eye(n, dtype=np.complex128)
    for lo, g in seq:                             # build the product in application order: M <- E M
        E = np.eye(n, dtype=np.complex128)
        E[lo:lo + 2, lo:lo + 2] = g
        Vd = E @ Vd
    Cfw = Vd[:, np.nonzero(occ)[0]]
    if Cfw.shape[1] != m:
        infid = 1.0
    else:
        s = np.linalg.svd(C.conj().T @ Cfw, compute_uv=False)
        infid = float(1.0 - np.prod(s) ** 2)
    return occ, seq, infid


def schedule_layers(seq):
    """ASAP layering of 2-site gates on bonds (lo, lo+1): returns list of layers, each a list of seq indices sorted
    by bond; gates in one layer act on disjoint sites, dependencies are preserved."""
    last = {}
    layer_of = []
    for k, (lo, _) in enumerate(seq):
        L = max(last.get(lo, -1), last.get(lo + 1, -1)) + 1
        layer_of.append(L)
        last[lo] = last[lo + 1] = L
    nL = max(layer_of) + 1 if layer_of else 0
    layers = [[] for _ in range(nL)]
    for k, L in enumerate(layer_of):
        layers[L].append(k)
    for L in layers:
        L.sort(key=lambda k: seq[k][0])
    return layers


EIG_STATS = {"fallback": 0}
SMALL_EIG_CPU = 48          # Hermitian eigenproblems up to this size go to the CPU (GPU eigh is latency-bound)


def _herm_eig(G: torch.Tensor):
    """Eigen-decomposition of a Hermitian PSD Gram matrix; GPU eigh for large blocks, CPU (complex128) for small
    blocks and as a fallback when cuSOLVER fails to converge (repeated zero eigenvalues in FP32)."""
    G = 0.5 * (G + G.mH)
    n = G.shape[0]
    if not bool(torch.isfinite(G).all()):
        raise FloatingPointError("non-finite Gram matrix")
    if G.device.type == "cuda" and n > SMALL_EIG_CPU:
        try:
            w, V = torch.linalg.eigh(G)
            if bool(torch.isfinite(w).all()) and bool(torch.isfinite(V).all()):
                return w, V
        except Exception:  # noqa: BLE001
            pass
        EIG_STATS["fallback"] += 1
    Gc = G.detach().to("cpu", torch.complex128).resolve_conj()
    w, V = torch.linalg.eigh(Gc)
    return w.to(G.device, G.real.dtype), V.to(G.device, G.dtype)


# --------------------------------------------------------------------------------------------------------- MPS

class SymMPS:
    """MPS with U(1)xU(1) bond labels; tensors A[i] (Dl, 4, Dr); single orthogonality centre."""

    def __init__(self, norb: int, nelec: tuple[int, int], device="cpu", dtype=torch.complex128,
                 occ: list[int] | None = None):
        self.n = norb
        self.nelec = tuple(nelec)
        self.device = torch.device(device)
        self.dtype = dtype
        self.site_q = torch.tensor(SITE_Q, dtype=torch.int64, device=self.device)
        if occ is None:   # Hartree-Fock product state in ffsim order: alpha 0..na-1, beta 0..nb-1 occupied
            occ = [int(p < nelec[0]) + 2 * int(p < nelec[1]) for p in range(norb)]
        self.A, self.q = [], [torch.zeros(1, dtype=torch.int64, device=self.device)]
        for p, s in enumerate(occ):
            t = torch.zeros((1, 4, 1), dtype=dtype, device=self.device)
            t[0, s, 0] = 1.0
            self.A.append(t)
            self.q.append(self.q[-1] + SITE_Q[s])
        self.center = 0
        self.discarded = 0.0           # accumulated discarded weight (sum over truncations)
        self.max_disc = 0.0
        self.n_trunc = 0
        self.stats = {"t_eigh": 0.0, "t_qr": 0.0, "n_gate": 0, "n_qr": 0}

    # ---------------------------------------------------------------- helpers
    def bond_dims(self):
        return [a.shape[2] for a in self.A[:-1]]

    def _row_q(self, i):          # labels of (l, s) rows of site i
        return (self.q[i][:, None] + self.site_q[None, :]).reshape(-1)

    def _col_q(self, i):          # labels of (s, r) cols of site i (left-flow label of the left bond)
        return (self.q[i + 1][None, :] - self.site_q[:, None]).reshape(-1)

    # ---------------------------------------------------------------- centre moves (block QR)
    def _qr_right(self, i):
        """site i -> left isometry, R absorbed into site i+1."""
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl * 4, Dr)
        qrow, qcol = self._row_q(i), self.q[i + 1]
        blocks_Q, blocks_R, labs, rows_l, cols_l = [], [], [], [], []
        for c in torch.unique(qcol).tolist():
            C = (qcol == c).nonzero().squeeze(1)
            Rr = (qrow == c).nonzero().squeeze(1)
            if len(Rr) == 0:
                continue
            Q, R = torch.linalg.qr(M[Rr][:, C])
            blocks_Q.append(Q)
            blocks_R.append(R)
            labs.append(c)
            rows_l.append(Rr)
            cols_l.append(C)
        K = sum(Q.shape[1] for Q in blocks_Q)
        Qf = torch.zeros((Dl * 4, K), dtype=self.dtype, device=self.device)
        Rf = torch.zeros((K, Dr), dtype=self.dtype, device=self.device)
        newq = torch.empty(K, dtype=torch.int64, device=self.device)
        o = 0
        for Q, R, c, Rr, C in zip(blocks_Q, blocks_R, labs, rows_l, cols_l):
            k = Q.shape[1]
            Qf[Rr, o:o + k] = Q
            Rf[o:o + k][:, C] = R
            newq[o:o + k] = c
            o += k
        self.A[i] = Qf.reshape(Dl, 4, K)
        self.A[i + 1] = torch.einsum("kr,rsb->ksb", Rf, self.A[i + 1])
        self.q[i + 1] = newq
        self.stats["t_qr"] += time.time() - t0
        self.stats["n_qr"] += 1

    def _qr_left(self, i):
        """site i -> right isometry, L absorbed into site i-1."""
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        qrow, qcol = self.q[i], self._col_q(i)
        blocks, o = [], 0
        for c in torch.unique(qrow).tolist():
            Rr = (qrow == c).nonzero().squeeze(1)
            C = (qcol == c).nonzero().squeeze(1)
            if len(C) == 0:
                continue
            Q, R = torch.linalg.qr(M[Rr][:, C].mH)          # M_c^dag = Q R  ->  M_c = R^dag Q^dag
            blocks.append((Q, R, c, Rr, C))
        K = sum(b[0].shape[1] for b in blocks)
        Qf = torch.zeros((4 * Dr, K), dtype=self.dtype, device=self.device)
        Lf = torch.zeros((Dl, K), dtype=self.dtype, device=self.device)
        newq = torch.empty(K, dtype=torch.int64, device=self.device)
        for Q, R, c, Rr, C in blocks:
            k = Q.shape[1]
            Qf[C, o:o + k] = Q
            Lf[Rr, o:o + k] = R.mH
            newq[o:o + k] = c
            o += k
        self.A[i] = Qf.mH.reshape(K, 4, Dr)
        self.A[i - 1] = torch.einsum("asl,lk->ask", self.A[i - 1], Lf)
        self.q[i] = newq
        self.stats["t_qr"] += time.time() - t0
        self.stats["n_qr"] += 1

    def move_center(self, target: int):
        while self.center < target:
            self._qr_right(self.center)
            self.center += 1
        while self.center > target:
            self._qr_left(self.center)
            self.center -= 1


    def _truncate(self, M, qrow, qcol, max_bond, cutoff, account: bool = True):
        """Block-wise (U(1)xU(1)) dominant left subspace of M (rows labelled qrow, cols qcol; M[r, c] = 0 unless
        qrow[r] == qcol[c]).  Returns (X isometry rows x k, labels (k,), kept weight, total weight)."""
        ws, vecs, labs, idxs = [], [], [], []
        common = set(torch.unique(qrow).tolist()) & set(torch.unique(qcol).tolist())
        for c in sorted(common):
            Rr = (qrow == c).nonzero().squeeze(1)
            C = (qcol == c).nonzero().squeeze(1)
            T = M[Rr][:, C]
            w, V = _herm_eig(T @ T.mH)
            ws.append(w.real)
            vecs.append(V)
            labs.append(c)
            idxs.append(Rr)
        allw = torch.cat(ws)
        total = float(allw.clamp(min=0).sum())
        order = torch.argsort(allw, descending=True)
        sw = allw[order]
        keep = int((sw > cutoff * total).sum()) if cutoff > 0 else len(sw)
        keep = max(1, min(keep, max_bond))
        kept_w = float(sw[:keep].clamp(min=0).sum())
        if account:
            disc = max(0.0, 1.0 - kept_w / total) if total > 0 else 0.0
            self.discarded += disc
            self.max_disc = max(self.max_disc, disc)
            self.n_trunc += 1
        sel_np = order[:keep].cpu().numpy()
        offs = np.cumsum([0] + [len(w) for w in ws])
        blk = np.searchsorted(offs, sel_np, side="right") - 1
        labs_np = np.array(labs)[blk]
        perm = np.lexsort((np.arange(keep), labs_np))
        sel_np, blk = sel_np[perm], blk[perm]
        newq = torch.as_tensor(labs_np[perm], dtype=torch.int64, device=self.device)
        X = torch.zeros((M.shape[0], keep), dtype=self.dtype, device=self.device)
        for b in np.unique(blk):
            cols = np.nonzero(blk == b)[0]
            loc = sel_np[cols] - offs[b]
            X[idxs[b][:, None], torch.as_tensor(cols, device=self.device)[None, :]] = \
                vecs[b][:, torch.as_tensor(loc, device=self.device)]
        return X, newq, kept_w, total

    def split_left(self, i: int, max_bond: int, cutoff: float = 0.0):
        """Centre at i -> truncate bond (i-1, i) optimally (canonical form), site i right isometry, centre -> i-1."""
        assert self.center == i and i > 0
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        Y, newq, kept_w, total = self._truncate(M.mH, self._col_q(i), self.q[i], max_bond, cutoff)
        k = Y.shape[1]
        self.A[i] = Y.mH.reshape(k, 4, Dr)
        self.A[i - 1] = torch.einsum("asl,lk->ask", self.A[i - 1], M @ Y)
        self.q[i] = newq
        self.center = i - 1

    def normalize_center(self):
        c = self.center
        self.A[c] = self.A[c] / torch.linalg.vector_norm(self.A[c])

    def apply_mpo_window(self, l: int, r: int, Ws, bl, br, dq, max_bond: int, cutoff: float = 0.0,
                         zip_margin: float = 2.0):
        """Exact windowed MPO (identity outside [l, r]) applied by zip-up (bond <= zip_margin * max_bond, no
        relative cutoff: the zip-up Gram is not a Schmidt spectrum), then compressed optimally by a canonical
        right-to-left sweep to max_bond / cutoff.  Ends with the centre at l."""
        self._zipup(l, r, Ws, bl, br, dq, int(zip_margin * max_bond))
        t0 = time.time()
        for i in range(r, l, -1):
            self.split_left(i, max_bond, cutoff)
        self.normalize_center()
        self.stats["t_compress"] = self.stats.get("t_compress", 0.0) + time.time() - t0

    def _zipup(self, l: int, r: int, Ws, bl, br, dq, max_bond: int, cutoff: float = 1e-15):
        """Zip-up application of an MPO that is the identity outside sites [l, r].
        Ws[i - l]: (D, 4, 4, D) tensors W[w, s_out, s_in, w']; bl, br: boundary vectors (D,); dq (D,) label shift
        (codes) carried by each channel (net (N_a, N_b) created by the operator part left of the bond).
        Ends with the centre at r."""
        t0 = time.time()
        self.move_center(l)
        D = len(bl)
        bl = torch.as_tensor(np.asarray(bl), dtype=self.dtype, device=self.device)
        br = torch.as_tensor(np.asarray(br), dtype=self.dtype, device=self.device)
        dq = torch.as_tensor(np.asarray(dq), dtype=torch.int64, device=self.device)
        chi_l = self.A[l].shape[0]
        C = torch.eye(chi_l, dtype=self.dtype, device=self.device)[:, None, :] * bl[None, :, None]
        qleft = self.q[l]
        for i in range(l, r + 1):
            W = torch.as_tensor(np.asarray(Ws[i - l]), dtype=self.dtype, device=self.device)
            A = self.A[i]
            CA = torch.einsum("xwa,asb->xwsb", C, A)
            T = torch.einsum("xwsb,wtsv->xtvb", CA, W)            # (chi', 4, D, chi_b)
            if i < r:
                x, _, _, cb = T.shape
                M = T.reshape(x * 4, D * cb)
                qrow = (qleft[:, None] + self.site_q[None, :]).reshape(-1)
                qcol = (dq[:, None] + self.q[i + 1][None, :]).reshape(-1)
                X, newq, kept_w, total = self._truncate(M, qrow, qcol, max_bond, cutoff, account=False)
                k = X.shape[1]
                self.A[i] = X.reshape(x, 4, k)
                C = (X.mH @ M).reshape(k, D, cb)
                self.q[i + 1] = newq
                qleft = newq
            else:
                T = torch.einsum("xtvb,v->xtb", T, br)
                nrm = torch.linalg.vector_norm(T)
                self.A[i] = T / nrm
        self.center = r
        self.stats["t_mpo"] = self.stats.get("t_mpo", 0.0) + time.time() - t0
        self.stats["n_mpo"] = self.stats.get("n_mpo", 0) + 1

    # ---------------------------------------------------------------- gates
    def apply_1site(self, i: int, diag4):
        d = torch.as_tensor(np.asarray(diag4), dtype=self.dtype, device=self.device)
        self.A[i] = self.A[i] * d[None, :, None]

    def apply_2site(self, i: int, G16, direction: str, max_bond: int, cutoff: float = 0.0, diag: bool = False):
        """Apply a 16x16 gate (or its diagonal if diag) on sites (i, i+1); the centre must be at i or i+1.
        direction 'right': site i becomes a left isometry, centre -> i+1; 'left': site i+1 right isometry, centre -> i."""
        assert self.center in (i, i + 1), (self.center, i)
        t0 = time.time()
        A, B = self.A[i], self.A[i + 1]
        Dl, Dr = A.shape[0], B.shape[2]
        th = torch.einsum("asm,mtb->astb", A, B)
        G = torch.as_tensor(np.asarray(G16), dtype=self.dtype, device=self.device)
        if diag:
            th = th * G.reshape(1, 4, 4, 1)
        else:
            th = torch.einsum("xy,ayb->axb", G, th.reshape(Dl, 16, Dr)).reshape(Dl, 4, 4, Dr)
        M = th.reshape(Dl * 4, 4 * Dr)
        qrow = self._row_q(i)
        qcol = self._col_q(i + 1)
        if direction == "right":
            X, newq, kept_w, total = self._truncate(M, qrow, qcol, max_bond, cutoff)
        else:
            Xc, newq, kept_w, total = self._truncate(M.mH, qcol, qrow, max_bond, cutoff)
            X = Xc
        keep = X.shape[1]
        nrm = math.sqrt(kept_w) if kept_w > 0 else 1.0
        if direction == "right":
            self.A[i] = X.reshape(Dl, 4, keep)
            self.A[i + 1] = (X.mH @ M).reshape(keep, 4, Dr) / nrm
            self.center = i + 1
        else:
            self.A[i + 1] = X.mH.reshape(keep, 4, Dr)
            self.A[i] = (M @ X).reshape(Dl, 4, keep) / nrm
            self.center = i
        self.q[i + 1] = newq
        self.stats["t_eigh"] += time.time() - t0
        self.stats["n_gate"] += 1

    # ---------------------------------------------------------------- circuit pieces
    def apply_gate_sequence(self, seq, max_bond, cutoff=0.0):
        """seq: [(lo, g 2x2)] in application order."""
        layers = schedule_layers(seq)
        right = self.center <= self.n // 2
        for L in layers:
            ks = L if right else L[::-1]
            for k in ks:
                lo, g = seq[k]
                if right:
                    self.move_center(lo)
                    self.apply_2site(lo, givens_gate16(g), "right", max_bond, cutoff)
                else:
                    self.move_center(lo + 1)
                    self.apply_2site(lo, givens_gate16(g), "left", max_bond, cutoff)
            right = not right

    def apply_phases(self, phases):
        for p, ph in enumerate(phases):
            self.apply_1site(p, [1.0, ph, ph, ph * ph])

    def apply_orbital_rotation(self, R, max_bond, cutoff=0.0):
        seq, phases = givens_sequence(R)
        self.apply_gate_sequence(seq, max_bond, cutoff)
        self.apply_phases(phases)

    def apply_diag_coulomb(self, z_nn, z_on, max_bond, cutoff=0.0):
        for p, z in enumerate(z_on):
            self.apply_1site(p, [1.0, 1.0, 1.0, np.exp(1j * z)])
        right = self.center <= self.n // 2
        bonds = range(self.n - 1) if right else range(self.n - 2, -1, -1)
        for p in bonds:
            if right:
                self.move_center(p)
                self.apply_2site(p, nn_dc_diag16(z_nn[p]), "right", max_bond, cutoff, diag=True)
            else:
                self.move_center(p + 1)
                self.apply_2site(p, nn_dc_diag16(z_nn[p]), "left", max_bond, cutoff, diag=True)

    # ---------------------------------------------------------------- dense (tiny systems, tests)
    def to_dense(self) -> np.ndarray:
        """Full 4^n vector, index = sum_p s_p 4^(n-1-p) (site 0 most significant)."""
        v = self.A[0].reshape(4, -1)
        for i in range(1, self.n):
            v = torch.einsum("xm,msb->xsb", v, self.A[i]).reshape(-1, self.A[i].shape[2])
        return v.reshape(-1).resolve_conj().cpu().numpy()


# ---------------------------------------------------------------------------- dense Hamiltonian (tiny systems)

def dense_fock_hamiltonian(one_body, two_body, constant, norb):
    """H in the MPS convention (JW order 0a,0b,1a,1b,...; index sum_p s_p 4^(n-1-p)), scipy sparse, 4^n x 4^n."""
    import scipy.sparse as sp
    nm = 2 * norb
    dim = 1 << nm
    # mode m = 2p + sigma; basis index bits: site p state s_p = n_pa + 2 n_pb at position 4^(n-1-p)

    def bit_of(m):
        p, sg = divmod(m, 2)
        return 2 * (norb - 1 - p) + sg

    idx = np.arange(dim)
    occ = np.array([(idx >> bit_of(m)) & 1 for m in range(nm)])        # (nm, dim)
    ops = []
    for m in range(nm):
        before = occ[:m].sum(0) if m > 0 else np.zeros(dim, dtype=int)
        ok = occ[m] == 0
        src = idx[ok]
        dst = src | (1 << bit_of(m))
        val = (-1.0) ** before[ok]
        ops.append(sp.csr_matrix((val, (dst, src)), shape=(dim, dim)))
    H = sp.csr_matrix((dim, dim), dtype=np.complex128)
    for p in range(norb):
        for q in range(norb):
            if abs(one_body[p, q]) < 1e-14:
                continue
            for s in range(2):
                H = H + one_body[p, q] * (ops[2 * p + s] @ ops[2 * q + s].T)
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for s_ in range(norb):
                    v = two_body[p, q, r, s_]
                    if abs(v) < 1e-14:
                        continue
                    for a in range(2):
                        for b in range(2):
                            H = H + 0.5 * v * (ops[2 * p + a] @ ops[2 * r + b] @ ops[2 * s_ + b].T
                                               @ ops[2 * q + a].T)
    return H + constant * sp.identity(dim, format="csr")


# ------------------------------------------------------------------------------------------- block2 expectation

class Block2Energy:
    """<mps|H|mps> with block2 (SZ, complex): the MPS is converted block by block (U(1)xU(1) labels -> SZ(n, 2Sz)),
    H is block2's quantum-chemistry MPO built from (complex) rotated integrals.  CPU, `n_threads` OpenMP threads."""

    def __init__(self, norb: int, nelec: tuple[int, int], scratch: str, n_threads: int = 8,
                 stack_mem: int = 4 << 30, sign_ab: float = 1.0):
        import block2  # noqa: F401
        from pyblock2.driver.core import DMRGDriver, SymmetryTypes
        self.norb, self.nelec = norb, tuple(nelec)
        self.driver = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SZ | SymmetryTypes.CPX,
                                 n_threads=n_threads, stack_mem=stack_mem)
        self.driver.initialize_system(n_sites=norb, n_elec=sum(nelec), spin=nelec[0] - nelec[1])
        self.sign_ab = sign_ab
        self._ntag = 0

    def mpo(self, h1e, g2e, ecore, algo_type=None):
        kw = {} if algo_type is None else {"algo_type": algo_type}
        return self.driver.get_qc_mpo(h1e=np.asarray(h1e, dtype=np.complex128),
                                      g2e=np.asarray(g2e, dtype=np.complex128), ecore=ecore, iprint=0, **kw)

    def to_block2(self, mps: SymMPS, tag: str | None = None):
        import block2 as b
        import block2.cpx as bx
        import block2.cpx.sz as bs
        import block2.sz as brs
        n = mps.n
        mps.move_center(0)
        if tag is None:
            tag = "TNKET"            # one tag, overwritten every call: no scratch growth over many reward calls

        def SZ(code):
            na, nb = divmod(int(code), QB)
            return b.SZ(na + nb, na - nb, 0)

        vacuum = b.SZ(0, 0, 0)
        target = b.SZ(sum(self.nelec), self.nelec[0] - self.nelec[1], 0)
        basis = []
        for _ in range(n):
            p = brs.StateInfo()
            p.allocate(4)
            for ix, c in enumerate(SITE_Q):
                p.quanta[ix] = SZ(c)
                p.n_states[ix] = 1
            p.sort_states()
            basis.append(p)
        info = brs.MPSInfo(n, vacuum, target, brs.VectorStateInfo(basis))
        info.tag = tag
        info.set_bond_dimension_full_fci(vacuum, vacuum)
        info.left_dims[0] = brs.StateInfo(vacuum)
        qs = [np.asarray(q) if isinstance(q, np.ndarray) else q.cpu().numpy() for q in mps.q]
        for bnd in range(1, n):
            labs, cnts = np.unique(qs[bnd], return_counts=True)
            p = info.left_dims[bnd]
            p.allocate(len(labs))
            for ix, (c, v) in enumerate(zip(labs, cnts)):
                p.quanta[ix] = SZ(c)
                p.n_states[ix] = int(v)
            p.sort_states()
            p = info.right_dims[bnd]
            p.allocate(len(labs))
            for ix, (c, v) in enumerate(zip(labs, cnts)):
                p.quanta[ix] = target - SZ(c)
                p.n_states[ix] = int(v)
            p.sort_states()
        info.left_dims[n] = brs.StateInfo(target)
        info.right_dims[0] = brs.StateInfo(target)
        info.right_dims[n] = brs.StateInfo(vacuum)
        info.bond_dim = info.get_max_bond_dimension()
        info.save_mutable()
        info.save_data("%s/%s-mps_info.bin" % (b.Global.frame.save_dir, tag))
        tensors = [bs.SparseTensor() for _ in range(n)]
        sgn = np.array([1.0, 1.0, 1.0, self.sign_ab])
        for i in range(n):
            Ai = mps.A[i]
            A = (np.asarray(Ai, dtype=np.complex128) if isinstance(Ai, np.ndarray) else
                 Ai.to(torch.complex128).resolve_conj().cpu().numpy()) * sgn[None, :, None]
            bb = basis[i]
            tensors[i].data = bs.VectorVectorPSSTensor([bs.VectorPSSTensor() for _ in range(bb.n)])
            ql_all, qr_all = qs[i], qs[i + 1]
            for s, cs in enumerate(SITE_Q):
                im = bb.find_state(SZ(cs))
                for cl in np.unique(ql_all):
                    rows = np.nonzero(ql_all == cl)[0]
                    cols = np.nonzero(qr_all == cl + cs)[0]
                    if len(cols) == 0:
                        continue
                    blk = np.ascontiguousarray(A[np.ix_(rows, [s], cols)])
                    if not np.any(blk):
                        continue
                    qlab = SZ(cl) if i > 0 else vacuum
                    qrab = SZ(cl + cs) if i < n - 1 else target
                    tensors[i].data[im].append(((qlab, qrab), bx.Tensor(b.VectorMKLInt(list(blk.shape)))))
                    np.array(tensors[i].data[im][-1][1], copy=False)[:] = blk
        umps = bs.UnfusedMPS()
        umps.info = info
        umps.n_sites = n
        umps.canonical_form = "K" + "R" * (n - 1)
        umps.center = 0
        umps.dot = 1
        umps.tensors = bs.VectorSpTensor(tensors)
        return umps.finalize()

    def expectation(self, bmps, mpo):
        return complex(self.driver.expectation(bmps, mpo, bmps))

    def norm2(self, bmps):
        return complex(self.driver.expectation(bmps, self.driver.get_identity_mpo(), bmps))


# ------------------------------------------------------------------------------------------------- public API

class LUCJEnergyTNGivens:
    """[Diagnostic, NOT usable as a reward] LUCJ energy by TEBD of the Givens-decomposed circuit in the ffsim (frame
    chain) orbital order.  The Clements networks pass through volume-law intermediate states (norb 15: discarded weight
    0.29 at chi 128 for one rotation), see the module docstring.  Kept for the record.

    LUCJ variational energy <psi|H|psi> from a truncated U(1)xU(1) MPS (TEBD on `device`) + block2 expectation.

    Mirrors pretrain.rl.energy.exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1)).
        ev = LUCJEnergyTNGivens(one_body, two_body, constant, norb, nelec, max_bond=512, device="cuda")
        E, info = ev.energy(U, Z, t1)
    info: discarded weight (sum over truncations, max single step), max bond, gate count, timings.
    """

    def __init__(self, one_body, two_body, constant, norb, nelec, max_bond: int = 256, cutoff: float = 1e-14,
                 device: str = "cuda", dtype=torch.complex64, block2_threads: int = 8, scratch: str | None = None,
                 stack_mem: int = 8 << 30, slater_init: bool = True, unitarize: bool = True):
        import tempfile
        self.h, self.g = np.asarray(one_body, dtype=np.float64), np.asarray(two_body, dtype=np.float64)
        self.const = float(constant)
        self.norb, self.nelec = int(norb), tuple(int(x) for x in nelec)
        self.max_bond, self.cutoff = int(max_bond), float(cutoff)
        self.device, self.dtype = device, dtype
        self.slater_init = slater_init and self.nelec[0] == self.nelec[1]
        self.unitarize = unitarize
        self._scratch = scratch or tempfile.mkdtemp(prefix="tn_b2_")
        self.b2 = Block2Energy(self.norb, self.nelec, self._scratch, n_threads=block2_threads, stack_mem=stack_mem)

    def state(self, U, Z, t1=None, max_bond: int | None = None):
        """Truncated MPS of exp(iJ_1) O(W) exp(iJ_0) O(U_0^dag)|HF> and the rotation F absorbed into H."""
        chi = self.max_bond if max_bond is None else int(max_bond)
        rots, dcs, F = lucj_layers(U, Z, t1, unitarize=self.unitarize)
        mps = SymMPS(self.norb, self.nelec, device=self.device, dtype=self.dtype)
        for k, (R, (znn, zon)) in enumerate(zip(rots, dcs)):
            if k == 0 and self.slater_init:
                mps.apply_gate_sequence(slater_sequence(R, self.nelec[0]), chi, self.cutoff)
            else:
                mps.apply_orbital_rotation(R, chi, self.cutoff)
            mps.apply_diag_coulomb(znn, zon, chi, self.cutoff)
        return mps, F

    def energy(self, U, Z, t1=None, max_bond: int | None = None):
        t0 = time.time()
        mps, F = self.state(U, Z, t1, max_bond)
        if mps.device.type == "cuda":
            torch.cuda.synchronize()
        t1_ = time.time()
        mps.move_center(0)
        nrm2 = float((mps.A[0].abs() ** 2).sum())
        h, g = rotate_hamiltonian(self.h, self.g, F)
        t2 = time.time()
        mpo = self.b2.mpo(h, g, self.const)
        t3 = time.time()
        bm = self.b2.to_block2(mps)
        t4 = time.time()
        e = self.b2.expectation(bm, mpo) / nrm2
        t5 = time.time()
        info = {"discarded_sum": mps.discarded, "discarded_max": mps.max_disc, "n_trunc": mps.n_trunc,
                "max_bond": max(mps.bond_dims()), "bond_dims": mps.bond_dims(), "imag": e.imag,
                "t_state": t1_ - t0, "t_rot_ints": t2 - t1_, "t_mpo": t3 - t2, "t_convert": t4 - t3,
                "t_expect": t5 - t4, "t_total": t5 - t0, **{k: v for k, v in mps.stats.items()}}
        del bm, mpo
        return float(e.real), info


# ------------------------------------------------------------------- split-localized basis: rotated-mode factors
#
# In an orbital basis S that does not mix occupied and virtual MOs (e.g. Boys-localized occupied + Boys-localized
# virtual orbitals, ordered along the molecule), |HF> is a product state and the LUCJ correlation is local, so the
# MPS stays small (norb 15: chi(1e-5 per cut) = 159 vs 1123 in the MO energy order).  The state in S coordinates is
#     phi_S = exp(i J_1(n^{V_1})) exp(i J_0(n^{V_0})) |HF_S>,   V_k = S^dag U_k,
# where n^{V}_p is the number operator of mode p = column p of V.  Every square-mask term is a commuting factor
#     exp(i z n_A n_B) = 1 + (e^{iz} - 1) b^dag_A b^dag_B b_B b_A        (A != B modes)
# applied exactly as a bond-dimension-17 MPO on the support window of the two modes (zip-up truncation).
# The final orbital rotation is absorbed into H: integrals rotated by final(t1) @ S.

_P4 = np.diag([1.0, -1.0, -1.0, 1.0])
_I4 = np.eye(4)
_CDAG_SITE = (np.array([[0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0]], dtype=float),     # c^dag_alpha
              np.array([[0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, -1, 0, 0]], dtype=float))    # c^dag_beta


def _mode_mpo(coef, spin, dagger):
    """MPO site tensors (L, 2, 4, 4, 2) of sum_q coef_q c^(dag)_{q,spin} with Jordan-Wigner strings."""
    op = _CDAG_SITE[spin] if dagger else _CDAG_SITE[spin].T
    W = np.zeros((len(coef), 2, 4, 4, 2), dtype=np.complex128)
    W[:, 0, :, :, 0] = _P4
    W[:, 0, :, :, 1] = np.asarray(coef)[:, None, None] * op[None]
    W[:, 1, :, :, 1] = _I4
    return W


def pair_factor_mpo(vA, sA, vB, sB, z, tol=1e-8):
    """exp(i z n_A n_B) for orthonormal modes A=(vA, spin sA), B=(vB, sB), A != B, as a windowed MPO.
    Returns (l, r, Ws (r-l+1, 17, 4, 4, 17), bl, br, dq) for SymMPS.apply_mpo_window."""
    amp = np.maximum(np.abs(vA), np.abs(vB))
    sup = np.nonzero(amp > tol * amp.max())[0]
    l, r = int(sup[0]), int(sup[-1])
    sl = slice(l, r + 1)
    W1 = _mode_mpo(vA[sl], sA, True)                  # b^dag_A
    W2 = _mode_mpo(vB[sl], sB, True)                  # b^dag_B
    W3 = _mode_mpo(np.conj(vB[sl]), sB, False)        # b_B
    W4 = _mode_mpo(np.conj(vA[sl]), sA, False)        # b_A
    X = np.einsum("iastA,ibtuB,icuvC,idvwD->iabcdswABCD", W1, W2, W3, W4, optimize=True)
    L = r - l + 1
    X = X.reshape(L, 16, 4, 4, 16)
    Ws = np.zeros((L, 17, 4, 4, 17), dtype=np.complex128)
    Ws[:, 0, :, :, 0] = _I4
    Ws[:, 1:, :, :, 1:] = X
    c = np.exp(1j * z) - 1.0
    bl = np.zeros(17, dtype=np.complex128)
    bl[0], bl[1] = 1.0, c                              # channel 1 = (0,0,0,0)
    br = np.zeros(17, dtype=np.complex128)
    br[0], br[16] = 1.0, 1.0                           # channel 16 = (1,1,1,1)
    e = (QB, 1)                                        # code of one alpha / one beta electron
    dq = np.zeros(17, dtype=np.int64)
    for ch in range(16):
        c1, c2, c3, c4 = (ch >> 3) & 1, (ch >> 2) & 1, (ch >> 1) & 1, ch & 1
        dq[1 + ch] = (c1 - c4) * e[sA] + (c2 - c3) * e[sB]
    return l, r, Ws, bl, br, dq


def square_factors(V: np.ndarray, Zk: np.ndarray):
    """(modeA, spinA, modeB, spinB, z) for the square-mask diagonal-Coulomb layer in rotated modes V[:, p]."""
    n = V.shape[0]
    Zs = 0.5 * (Zk + Zk.T)
    out = []
    for p in range(n):
        if abs(Zs[p, p]) > 0:
            out.append((V[:, p], 0, V[:, p], 1, float(Zs[p, p])))
    for p in range(n - 1):
        z = float(Zs[p, p + 1])
        if abs(z) > 0:
            out.append((V[:, p], 0, V[:, p + 1], 0, z))
            out.append((V[:, p], 1, V[:, p + 1], 1, z))
    return out


class LUCJEnergySplitTN:
    """LUCJ energy with the MPS in a fixed orbital basis S (MO coordinates, occupied/virtual not mixed, ordered).

        ev = LUCJEnergySplitTN(one_body, two_body, constant, norb, nelec, S, occ_mask, max_bond=256, device="cuda")
        E, info = ev.energy(U, Z, t1)
    """

    def __init__(self, one_body, two_body, constant, norb, nelec, S, occ_mask, max_bond: int = 256,
                 cutoff: float = 1e-12, mode_tol: float = 1e-8, device: str = "cuda", dtype=torch.complex64,
                 block2_threads: int = 8, scratch: str | None = None, stack_mem: int = 8 << 30,
                 unitarize: bool = True):
        import tempfile
        assert nelec[0] == nelec[1], "closed-shell only"
        self.h, self.g = np.asarray(one_body, dtype=np.float64), np.asarray(two_body, dtype=np.float64)
        self.const = float(constant)
        self.norb, self.nelec = int(norb), tuple(int(x) for x in nelec)
        self.S = np.asarray(S)
        self.occ = np.asarray(occ_mask, dtype=bool)
        assert self.occ.sum() == self.nelec[0]
        self.max_bond, self.cutoff, self.mode_tol = int(max_bond), float(cutoff), float(mode_tol)
        self.device, self.dtype, self.unitarize = device, dtype, unitarize
        self._scratch = scratch or tempfile.mkdtemp(prefix="tn_b2_")
        self.b2 = Block2Energy(self.norb, self.nelec, self._scratch, n_threads=block2_threads, stack_mem=stack_mem)
        self._mpo_cache = {}
        self.zip_margin = 2.0

    def state(self, U, Z, t1=None, max_bond: int | None = None):
        chi = self.max_bond if max_bond is None else int(max_bond)
        U = np.asarray(U, dtype=np.complex128)
        if self.unitarize:
            U = np.stack([polar_unitary(u) for u in U])
        Z = np.asarray(Z, dtype=np.float64)
        occs = [3 if o else 0 for o in self.occ]
        if str(self.device) == "cpu":
            mps = NpSymMPS(self.norb, self.nelec, occs)
        else:
            mps = GpuSymMPS(self.norb, self.nelec, occs, device=self.device, dtype=self.dtype)
        for k in range(U.shape[0]):
            V = self.S.conj().T @ U[k]
            facs = []
            for vA, sA, vB, sB, z in square_factors(V, Z[k]):
                facs.append(pair_factor_mpo(vA, sA, vB, sB, z, self.mode_tol))
            facs.sort(key=lambda f: (f[0], f[1]))
            for l, r, Ws, bl, br, dq in facs:
                mps.apply_mpo_window(l, r, Ws, bl, br, dq, chi, self.cutoff, zip_margin=self.zip_margin)
        if t1 is not None:
            from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
            Fr = orbital_rotation_from_t1_amplitudes(np.asarray(t1, dtype=np.float64)) @ self.S
        else:
            Fr = self.S
        return mps, Fr

    def energy(self, U, Z, t1=None, max_bond: int | None = None):
        t0 = time.time()
        mps, Fr = self.state(U, Z, t1, max_bond)
        if mps.device.type == "cuda":
            torch.cuda.synchronize()
        t1_ = time.time()
        mps.move_center(0)
        nrm2 = float((abs(mps.A[0]) ** 2).sum())
        key = hash(np.asarray(Fr).tobytes())  # Fr = final(t1) @ S: one MPO per molecule
        if key not in self._mpo_cache:
            self._mpo_cache.clear()
            h, g = rotate_hamiltonian(self.h, self.g, Fr)
            self._mpo_cache[key] = self.b2.mpo(h, g, self.const)
        mpo = self._mpo_cache[key]
        t2 = time.time()
        bm = self.b2.to_block2(mps)
        t3 = time.time()
        e = self.b2.expectation(bm, mpo) / nrm2
        t4 = time.time()
        info = {"discarded_sum": mps.discarded, "discarded_max": mps.max_disc, "n_trunc": mps.n_trunc,
                "max_bond": max(mps.bond_dims()), "imag": e.imag, "t_state": t1_ - t0, "t_mpo_build": t2 - t1_,
                "t_convert": t3 - t2, "t_expect": t4 - t3, "t_total": t4 - t0,
                **{k: v for k, v in mps.stats.items()}}
        del bm
        return float(e.real), info


def random_split_basis(n, nocc, rng):
    """Random occupied/virtual-preserving orthogonal basis with a random site order (tests)."""
    def rorth(m):
        Q, R = np.linalg.qr(rng.normal(size=(m, m)))
        return Q * np.sign(np.diagonal(R))
    S = np.zeros((n, n))
    S[:nocc, :nocc] = rorth(nocc)
    S[nocc:, nocc:] = rorth(n - nocc)
    perm = rng.permutation(n)
    occ = np.arange(n) < nocc
    return S[:, perm], occ[perm]


def split_localized_basis(name: str, root=None, method: str = "boys", order: str = "fiedler",
                          cache_dir: str | None = None):
    """Occupied/virtual split-localized orbitals of the active space, in MO coordinates (columns), ordered in 1D.

    method: 'boys' | 'pm' (pyscf.lo, AO geometry from jobs/<name>/<name>.xyz, active MOs from rhf_dataset).
    order : 'axis'    -- orbital centroids sorted along the principal axis of the centroid cloud;
            'fiedler' -- Fiedler vector of the graph with weights exp(-|R_i - R_j| / Angstrom) between centroids.
    Returns (S (n, n) real orthogonal, occ_mask (n,) bool).  The occupied/virtual blocks are re-orthonormalized
    (polar) so S is exactly block-diagonal in the MO basis of the Hamiltonian (|HF> is a product state)."""
    from pathlib import Path
    root = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    if cache_dir is not None:
        f = Path(cache_dir) / f"{name}_{method}_{order}.npz"
        if f.exists():
            d = np.load(f)
            return d["S"], d["occ"]
    from pyscf import gto, lo
    from pretrain.rl.hamiltonian import read_xyz
    mol = gto.M(atom=read_xyz(root / "jobs" / name / f"{name}.xyz"), basis="sto-3g", verbose=0)
    d = np.load(root / "rhf_dataset" / f"{name}.npz")
    C = d["mo_coeff"].astype(np.float64)
    nocc = int(d["nocc"])
    n = C.shape[1]
    Sao = mol.intor("int1e_ovlp")
    Loc = lo.Boys if method == "boys" else lo.PM
    S = np.zeros((n, n))
    for sl in (slice(0, nocc), slice(nocc, n)):
        L = Loc(mol, C[:, sl]).kernel()
        B = C[:, sl].T @ Sao @ L                       # (block, block), ~orthogonal (float32 MOs)
        S[sl, sl] = polar_unitary(B).real
    occ = np.arange(n) < nocc
    r = mol.intor("int1e_r")
    Lao = C @ S
    cen = np.einsum("xmn,mi,ni->ix", r, Lao, Lao) * 0.529177210903     # Bohr -> Angstrom
    if order == "axis":
        X = cen - cen.mean(0)
        w, V = np.linalg.eigh(X.T @ X)
        perm = np.argsort(X @ V[:, -1], kind="stable")
    elif order == "fiedler":
        D = np.linalg.norm(cen[:, None] - cen[None], axis=-1)
        Wt = np.exp(-D)
        np.fill_diagonal(Wt, 0.0)
        Lap = np.diag(Wt.sum(1)) - Wt
        ew, ev = np.linalg.eigh(Lap)
        perm = np.argsort(ev[:, 1], kind="stable")
    else:
        raise ValueError(order)
    S, occ = S[:, perm], occ[perm]
    if cache_dir is not None:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        np.savez(Path(cache_dir) / f"{name}_{method}_{order}.npz", S=S, occ=occ, centroids=cen[perm])
    return S, occ


# ------------------------------------------------------------------------------ numpy engine (CPU, low overhead)

SITE_Q_NP = np.array(SITE_Q, dtype=np.int64)


def _sectors(sorted_labels):
    """label -> (start, stop) for a sorted label array."""
    if len(sorted_labels) == 0:
        return {}
    cut = np.nonzero(np.diff(sorted_labels))[0] + 1
    starts = np.r_[0, cut]
    stops = np.r_[cut, len(sorted_labels)]
    return {int(sorted_labels[a]): (int(a), int(b)) for a, b in zip(starts, stops)}


class NpSymMPS:
    """numpy twin of SymMPS (complex128, CPU): U(1)xU(1) bond labels, single centre; sector bookkeeping by sorting so
    every block operation is a contiguous slice.  Used for small/moderate bond dimensions where GPU launch latency
    dominates (cf. pretrain/rl/tests profiling: ~14 sector eigh per truncation)."""

    def __init__(self, norb, nelec, occ, dtype=np.complex128):
        self.n, self.nelec, self.dtype = norb, tuple(nelec), dtype
        self.A, self.q = [], [np.zeros(1, dtype=np.int64)]
        for s in occ:
            t = np.zeros((1, 4, 1), dtype=dtype)
            t[0, s, 0] = 1.0
            self.A.append(t)
            self.q.append(self.q[-1] + SITE_Q_NP[s])
        self.center = 0
        self.discarded, self.max_disc, self.n_trunc = 0.0, 0.0, 0
        self.stats = {"t_zip": 0.0, "t_compress": 0.0, "t_qr": 0.0, "n_mpo": 0}
        self.zip_cutoff = 1e-15
        self.device = torch.device("cpu")

    def bond_dims(self):
        return [a.shape[2] for a in self.A[:-1]]

    # ------------------------------------------------------------ truncation core
    def _trunc(self, M, qrow, qcol, max_bond, cutoff, account=True):
        """Dominant row-space isometry of block-diagonal M (rows qrow, cols qcol). Returns X, labels, kept_w, total."""
        pr = np.argsort(qrow, kind="stable")
        pc = np.argsort(qcol, kind="stable")
        sr, sc = qrow[pr], qcol[pc]
        secr, secc = _sectors(sr), _sectors(sc)
        Ms = M[pr][:, pc]
        blocks = []
        for c, (r0, r1) in secr.items():
            if c not in secc:
                continue
            c0, c1 = secc[c]
            B = Ms[r0:r1, c0:c1]
            blocks.append((c, r0, B @ B.conj().T))
        ws, vecs, labs, offs = [None] * len(blocks), [None] * len(blocks), [], []
        for j, (c, r0, G) in enumerate(blocks):
            labs.append(c)
            offs.append(r0)
        # batch eigh: blocks of size <= 16 padded to the bin size, larger blocks individually
        small = [j for j, b in enumerate(blocks) if b[2].shape[0] <= 16]
        for binsz in (2, 4, 8, 16):
            idx = [j for j in small if binsz // 2 < blocks[j][2].shape[0] <= binsz or (binsz == 2 and blocks[j][2].shape[0] <= 2)]
            if not idx:
                continue
            stack = np.zeros((len(idx), binsz, binsz), dtype=M.dtype)
            for t, j in enumerate(idx):
                m = blocks[j][2].shape[0]
                stack[t, :m, :m] = blocks[j][2]
                if m < binsz:                                   # padded diagonal pushed below every real weight
                    stack[t, np.arange(m, binsz), np.arange(m, binsz)] = -1.0
            w_all, V_all = np.linalg.eigh(stack)
            for t, j in enumerate(idx):
                m = blocks[j][2].shape[0]
                keep_ = np.nonzero(np.abs(V_all[t, m:, :]).sum(0) < 1e-12)[0] if m < binsz else np.arange(binsz)
                keep_ = keep_[-m:] if len(keep_) >= m else keep_
                ws[j] = w_all[t][keep_]
                vecs[j] = V_all[t][:m, keep_]
        for j, (c, r0, G) in enumerate(blocks):
            if ws[j] is None:
                w, V = np.linalg.eigh(G)
                ws[j], vecs[j] = w, V
        allw = np.concatenate(ws)
        total = float(np.clip(allw, 0, None).sum())
        order = np.argsort(-allw, kind="stable")
        sw = allw[order]
        keep = int((sw > cutoff * total).sum()) if cutoff > 0 else len(sw)
        keep = max(1, min(keep, max_bond))
        kept_w = float(np.clip(sw[:keep], 0, None).sum())
        if account:
            disc = max(0.0, 1.0 - kept_w / total) if total > 0 else 0.0
            self.discarded += disc
            self.max_disc = max(self.max_disc, disc)
            self.n_trunc += 1
        bounds = np.cumsum([0] + [len(w) for w in ws])
        sel = order[:keep]
        blk = np.searchsorted(bounds, sel, side="right") - 1
        perm = np.lexsort((sel, blk))                  # group by block (block order = label order), then index
        sel, blk = sel[perm], blk[perm]
        X = np.zeros((M.shape[0], keep), dtype=M.dtype)
        labels = np.empty(keep, dtype=np.int64)
        j = 0
        for b in np.unique(blk):
            idx = sel[blk == b] - bounds[b]
            m = len(idx)
            r0 = offs[b]
            rows = pr[r0:r0 + vecs[b].shape[0]]
            X[rows, j:j + m] = vecs[b][:, idx]
            labels[j:j + m] = labs[b]
            j += m
        return X, labels, kept_w, total

    def _qr_right(self, i):
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl * 4, Dr)
        qrow = (self.q[i][:, None] + SITE_Q_NP[None, :]).reshape(-1)
        qcol = self.q[i + 1]
        pr = np.argsort(qrow, kind="stable")
        secr = _sectors(qrow[pr])
        secc = _sectors(qcol)                               # bond labels are kept sorted
        Qs, Rs, labs, rws, cls = [], [], [], [], []
        for c, (c0, c1) in secc.items():
            if c not in secr:
                continue
            r0, r1 = secr[c]
            rows = pr[r0:r1]
            Q, R = np.linalg.qr(M[rows, c0:c1])
            Qs.append(Q); Rs.append(R); labs.append(c); rws.append(rows); cls.append((c0, c1))
        K = sum(Q.shape[1] for Q in Qs)
        Qf = np.zeros((Dl * 4, K), dtype=A.dtype)
        Rf = np.zeros((K, Dr), dtype=A.dtype)
        newq = np.empty(K, dtype=np.int64)
        o = 0
        for Q, R, c, rows, (c0, c1) in zip(Qs, Rs, labs, rws, cls):
            k = Q.shape[1]
            Qf[rows, o:o + k] = Q
            Rf[o:o + k, c0:c1] = R
            newq[o:o + k] = c
            o += k
        self.A[i] = Qf.reshape(Dl, 4, K)
        B = self.A[i + 1]
        self.A[i + 1] = (Rf @ B.reshape(B.shape[0], -1)).reshape(K, 4, B.shape[2])
        self.q[i + 1] = newq
        self.stats["t_qr"] += time.time() - t0

    def _qr_left(self, i):
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        qrow = self.q[i]
        qcol = (self.q[i + 1][None, :] - SITE_Q_NP[:, None]).reshape(-1)
        pc = np.argsort(qcol, kind="stable")
        secc = _sectors(qcol[pc])
        secr = _sectors(qrow)
        parts, o = [], 0
        for c, (r0, r1) in secr.items():
            if c not in secc:
                continue
            c0, c1 = secc[c]
            cols = pc[c0:c1]
            Q, R = np.linalg.qr(M[r0:r1][:, cols].conj().T)
            parts.append((Q, R, c, (r0, r1), cols))
        K = sum(p[0].shape[1] for p in parts)
        Qf = np.zeros((4 * Dr, K), dtype=A.dtype)
        Lf = np.zeros((Dl, K), dtype=A.dtype)
        newq = np.empty(K, dtype=np.int64)
        for Q, R, c, (r0, r1), cols in parts:
            k = Q.shape[1]
            Qf[cols, o:o + k] = Q
            Lf[r0:r1, o:o + k] = R.conj().T
            newq[o:o + k] = c
            o += k
        self.A[i] = Qf.conj().T.reshape(K, 4, Dr)
        P = self.A[i - 1]
        self.A[i - 1] = (P.reshape(-1, Dl) @ Lf).reshape(P.shape[0], 4, K)
        self.q[i] = newq
        self.stats["t_qr"] += time.time() - t0

    def move_center(self, target):
        while self.center < target:
            self._qr_right(self.center)
            self.center += 1
        while self.center > target:
            self._qr_left(self.center)
            self.center -= 1

    def split_left(self, i, max_bond, cutoff=0.0):
        """Centre at i: optimal truncation of bond (i-1, i) via the (smaller) row Gram, centre -> i-1."""
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        qcol = (self.q[i + 1][None, :] - SITE_Q_NP[:, None]).reshape(-1)
        Ul, labels, kept_w, total = self._trunc(M, self.q[i], qcol, max_bond, cutoff)
        # M = U s Y^dag  ->  Y^dag = s^-1 U^dag M ; centre factor (absorbed left) = U s
        UM = Ul.conj().T @ M                                  # (k, 4Dr) = s Y^dag
        s = np.sqrt(np.clip(np.einsum("kc,kc->k", UM, UM.conj()).real, 1e-300, None))
        Yd = UM / s[:, None]
        k = Ul.shape[1]
        self.A[i] = Yd.reshape(k, 4, Dr)
        P = self.A[i - 1]
        self.A[i - 1] = (P.reshape(-1, Dl) @ (Ul * s[None, :])).reshape(P.shape[0], 4, k)
        self.q[i] = labels
        self.center = i - 1

    def normalize_center(self):
        c = self.center
        self.A[c] = self.A[c] / np.linalg.norm(self.A[c])

    def apply_mpo_window(self, l, r, Ws, bl, br, dq, max_bond, cutoff=0.0, zip_margin=2.0):
        t0 = time.time()
        self.move_center(l)
        D = len(bl)
        chi_l = self.A[l].shape[0]
        C = np.eye(chi_l, dtype=self.dtype)[:, None, :] * np.asarray(bl)[None, :, None]
        qleft = self.q[l]
        zmax = int(zip_margin * max_bond)
        for i in range(l, r + 1):
            W = Ws[i - l]
            A = self.A[i]
            x, _, a = C.shape
            b = A.shape[2]
            CA = (C.reshape(x * D, a) @ A.reshape(a, 4 * b)).reshape(x, D, 4, b)
            T = np.tensordot(CA, W, axes=([1, 2], [0, 2]))         # (x, b, t, D')
            if i < r:
                M = T.transpose(0, 2, 3, 1).reshape(x * 4, D * b)
                qrow = (qleft[:, None] + SITE_Q_NP[None, :]).reshape(-1)
                qcol = (np.asarray(dq)[:, None] + self.q[i + 1][None, :]).reshape(-1)
                X, newq, _, _ = self._trunc(M, qrow, qcol, zmax, self.zip_cutoff, account=False)
                k = X.shape[1]
                self.A[i] = X.reshape(x, 4, k)
                C = (X.conj().T @ M).reshape(k, D, b)
                self.q[i + 1] = newq
                qleft = newq
            else:
                Tf = np.tensordot(T, np.asarray(br), axes=([3], [0]))  # (x, b, t)
                self.A[i] = Tf.transpose(0, 2, 1) / np.linalg.norm(Tf)
        self.center = r
        t1 = time.time()
        for i in range(r, l, -1):
            self.split_left(i, max_bond, cutoff)
        self.normalize_center()
        self.stats["t_zip"] += t1 - t0
        self.stats["t_compress"] += time.time() - t1
        self.stats["n_mpo"] += 1

    def to_dense(self):
        v = self.A[0].reshape(4, -1)
        for i in range(1, self.n):
            v = (v @ self.A[i].reshape(self.A[i].shape[0], -1)).reshape(-1, self.A[i].shape[2])
        return v.reshape(-1)


# ------------------------------------------------------------------- GPU engine (labels on host, tensors on device)

def _batched_eigh_np(mats):
    """eigh of a list of small Hermitian numpy matrices; equal sizes are stacked into one LAPACK batch call."""
    out = [None] * len(mats)
    by = {}
    for j, G in enumerate(mats):
        by.setdefault(G.shape[0], []).append(j)
    for m, idx in by.items():
        if len(idx) == 1:
            out[idx[0]] = np.linalg.eigh(mats[idx[0]])
        else:
            w, V = np.linalg.eigh(np.stack([mats[j] for j in idx]))
            for t, j in enumerate(idx):
                out[j] = (w[t], V[t])
    return out


class GpuSymMPS(NpSymMPS):
    """SymMPS variant for large bond dimensions: tensors (complex64/128) on a CUDA device, U(1)xU(1) labels as numpy
    arrays on the host (no device syncs for bookkeeping).  Gram blocks are formed on the GPU; blocks up to
    `gpu_eig_min` are diagonalized on the CPU in one transfer (batched LAPACK), larger ones with GPU eigh."""

    gpu_eig_min = 256

    def __init__(self, norb, nelec, occ, device="cuda", dtype=torch.complex64):
        self.n, self.nelec = norb, tuple(nelec)
        self.device = torch.device(device)
        self.tdtype = dtype
        self.dtype = np.complex64 if dtype == torch.complex64 else np.complex128
        self.A, self.q = [], [np.zeros(1, dtype=np.int64)]
        for s in occ:
            t = torch.zeros((1, 4, 1), dtype=dtype, device=self.device)
            t[0, s, 0] = 1.0
            self.A.append(t)
            self.q.append(self.q[-1] + SITE_Q_NP[s])
        self.center = 0
        self.discarded, self.max_disc, self.n_trunc = 0.0, 0.0, 0
        self.stats = {"t_zip": 0.0, "t_compress": 0.0, "t_qr": 0.0, "n_mpo": 0}
        self.zip_cutoff = 1e-15

    def _t(self, x):
        return torch.as_tensor(x, device=self.device)

    def _trunc(self, M, qrow, qcol, max_bond, cutoff, account=True):
        pr = np.argsort(qrow, kind="stable")
        pc = np.argsort(qcol, kind="stable")
        secr, secc = _sectors(qrow[pr]), _sectors(qcol[pc])
        pr_t, pc_t = self._t(pr), self._t(pc)
        Ms = M.index_select(0, pr_t).index_select(1, pc_t)
        labs, offs, Gs = [], [], []
        for c, (r0, r1) in secr.items():
            if c not in secc:
                continue
            c0, c1 = secc[c]
            B = Ms[r0:r1, c0:c1]
            Gs.append(B @ B.mH)
            labs.append(c)
            offs.append(r0)
        small = [j for j, G in enumerate(Gs) if G.shape[0] < self.gpu_eig_min]
        res = [None] * len(Gs)
        if small:
            flat = torch.cat([Gs[j].reshape(-1) for j in small]).to(torch.complex128).cpu().numpy()
            mats, o = [], 0
            for j in small:
                m = Gs[j].shape[0]
                G = flat[o:o + m * m].reshape(m, m)
                mats.append(0.5 * (G + G.conj().T))
                o += m * m
            for j, wv in zip(small, _batched_eigh_np(mats)):
                res[j] = wv
        for j, G in enumerate(Gs):
            if res[j] is None:
                w, V = _herm_eig(G)
                res[j] = (w.double().cpu().numpy(), V)
        ws = [r[0] for r in res]
        allw = np.concatenate(ws)
        total = float(np.clip(allw, 0, None).sum())
        order = np.argsort(-allw, kind="stable")
        sw = allw[order]
        keep = int((sw > cutoff * total).sum()) if cutoff > 0 else len(sw)
        keep = max(1, min(keep, max_bond))
        kept_w = float(np.clip(sw[:keep], 0, None).sum())
        if account:
            disc = max(0.0, 1.0 - kept_w / total) if total > 0 else 0.0
            self.discarded += disc
            self.max_disc = max(self.max_disc, disc)
            self.n_trunc += 1
        bounds = np.cumsum([0] + [len(w) for w in ws])
        sel = order[:keep]
        blk = np.searchsorted(bounds, sel, side="right") - 1
        perm = np.lexsort((sel, blk))
        sel, blk = sel[perm], blk[perm]
        X = torch.zeros((M.shape[0], keep), dtype=M.dtype, device=self.device)
        labels = np.empty(keep, dtype=np.int64)
        j = 0
        host_parts, host_rows, host_cols = [], [], []
        for b in np.unique(blk):
            idx = sel[blk == b] - bounds[b]
            m = len(idx)
            V = res[b][1]
            rows = pr[offs[b]:offs[b] + V.shape[0]]
            if isinstance(V, np.ndarray):
                host_parts.append(V[:, idx])
                host_rows.append(rows)
                host_cols.append(np.arange(j, j + m))
            else:
                X[self._t(rows)[:, None], self._t(np.arange(j, j + m))[None, :]] = V[:, self._t(idx)].to(M.dtype)
            labels[j:j + m] = labs[b]
            j += m
        if host_parts:
            # one transfer: scatter all host-diagonalized blocks
            rr = np.concatenate([np.repeat(r, len(c)) for r, c in zip(host_rows, host_cols)])
            cc = np.concatenate([np.tile(c, len(r)) for r, c in zip(host_rows, host_cols)])
            vv = np.concatenate([p.reshape(-1) for p in host_parts])
            X[self._t(rr), self._t(cc)] = torch.as_tensor(vv, device=self.device).to(M.dtype)
        return X, labels, kept_w, total

    def _qr_right(self, i):
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl * 4, Dr)
        qrow = (self.q[i][:, None] + SITE_Q_NP[None, :]).reshape(-1)
        X, labels, _, _ = self._trunc(M, qrow, self.q[i + 1], 10 ** 9, 1e-14, account=False)
        k = X.shape[1]
        self.A[i] = X.reshape(Dl, 4, k)
        B = self.A[i + 1]
        self.A[i + 1] = (X.mH @ M @ B.reshape(Dr, -1)).reshape(k, 4, B.shape[2])
        self.q[i + 1] = labels
        self.stats["t_qr"] += time.time() - t0

    def _qr_left(self, i):
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        qcol = (self.q[i + 1][None, :] - SITE_Q_NP[:, None]).reshape(-1)
        Y, labels, _, _ = self._trunc(M.mH.contiguous(), qcol, self.q[i], 10 ** 9, 1e-14, account=False)
        k = Y.shape[1]
        self.A[i] = Y.mH.reshape(k, 4, Dr)
        P = self.A[i - 1]
        self.A[i - 1] = (P.reshape(-1, Dl) @ (M @ Y)).reshape(P.shape[0], 4, k)
        self.q[i] = labels
        self.stats["t_qr"] += time.time() - t0

    def split_left(self, i, max_bond, cutoff=0.0):
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        qcol = (self.q[i + 1][None, :] - SITE_Q_NP[:, None]).reshape(-1)
        Ul, labels, kept_w, total = self._trunc(M, self.q[i], qcol, max_bond, cutoff)
        UM = Ul.mH @ M
        s = torch.sqrt(torch.clamp((UM.abs() ** 2).sum(1), min=1e-30))
        Yd = UM / s[:, None]
        k = Ul.shape[1]
        self.A[i] = Yd.reshape(k, 4, Dr)
        P = self.A[i - 1]
        self.A[i - 1] = (P.reshape(-1, Dl) @ (Ul * s[None, :])).reshape(P.shape[0], 4, k)
        self.q[i] = labels
        self.center = i - 1

    def normalize_center(self):
        c = self.center
        self.A[c] = self.A[c] / torch.linalg.vector_norm(self.A[c])

    def apply_mpo_window(self, l, r, Ws, bl, br, dq, max_bond, cutoff=0.0, zip_margin=2.0):
        t0 = time.time()
        self.move_center(l)
        D = len(bl)
        chi_l = self.A[l].shape[0]
        blt = torch.as_tensor(np.asarray(bl), dtype=self.tdtype, device=self.device)
        brt = torch.as_tensor(np.asarray(br), dtype=self.tdtype, device=self.device)
        Wt = torch.as_tensor(np.asarray(Ws), dtype=self.tdtype, device=self.device)
        C = torch.eye(chi_l, dtype=self.tdtype, device=self.device)[:, None, :] * blt[None, :, None]
        qleft = self.q[l]
        zmax = int(zip_margin * max_bond)
        for i in range(l, r + 1):
            W = Wt[i - l]
            A = self.A[i]
            x, _, a = C.shape
            b = A.shape[2]
            CA = (C.reshape(x * D, a) @ A.reshape(a, 4 * b)).reshape(x, D, 4, b)
            T = torch.tensordot(CA, W, dims=([1, 2], [0, 2]))      # (x, b, t, D')
            if i < r:
                M = T.permute(0, 2, 3, 1).reshape(x * 4, D * b)
                qrow = (qleft[:, None] + SITE_Q_NP[None, :]).reshape(-1)
                qcol = (np.asarray(dq)[:, None] + self.q[i + 1][None, :]).reshape(-1)
                X, newq, _, _ = self._trunc(M, qrow, qcol, zmax, self.zip_cutoff, account=False)
                k = X.shape[1]
                self.A[i] = X.reshape(x, 4, k)
                C = (X.mH @ M).reshape(k, D, b)
                self.q[i + 1] = newq
                qleft = newq
            else:
                Tf = torch.tensordot(T, brt, dims=([3], [0]))       # (x, b, t)
                Tf = Tf.permute(0, 2, 1)
                self.A[i] = (Tf / torch.linalg.vector_norm(Tf)).contiguous()
        self.center = r
        t1 = time.time()
        for i in range(r, l, -1):
            self.split_left(i, max_bond, cutoff)
        self.normalize_center()
        if self.device.type == "cuda":
            torch.cuda.synchronize()
        self.stats["t_zip"] += t1 - t0
        self.stats["t_compress"] += time.time() - t1
        self.stats["n_mpo"] += 1

    def to_dense(self):
        v = self.A[0].reshape(4, -1)
        for i in range(1, self.n):
            v = (v @ self.A[i].reshape(self.A[i].shape[0], -1)).reshape(-1, self.A[i].shape[2])
        return v.reshape(-1).resolve_conj().cpu().numpy()


# --------------------------------------------------------------- integral-only split localization (ER) + ordering

def _er_localize_block(g: np.ndarray, idx: np.ndarray, sweeps: int = 50, tol: float = 1e-9):
    """Edmiston-Ruedenberg localization (maximize sum_k (kk|kk)) of the MOs `idx`, using only the MO-basis two-body
    tensor g (chemists' notation, real, 8-fold symmetric).  Jacobi sweeps with the exact 2x2 angle
        D(t) = const + A cos 4t + B sin 4t,  A = [(ii|ii)+(jj|jj)]/4 - [(ii|jj)+2(ij|ij)]/2,  B = (ii|ij)-(jj|ij),
    and O(m^3) slice updates of the integrals.  Returns R (m, m): columns = localized orbitals (block MO coords)."""
    m = len(idx)
    R = np.eye(m)
    if m < 2:
        return R
    gb = np.ascontiguousarray(g[np.ix_(idx, idx, idx, idx)])
    for _ in range(sweeps):
        delta = 0.0
        for i in range(m - 1):
            for j in range(i + 1, m):
                A = 0.25 * (gb[i, i, i, i] + gb[j, j, j, j]) - 0.5 * (gb[i, i, j, j] + 2.0 * gb[i, j, i, j])
                B = gb[i, i, i, j] - gb[j, j, i, j]
                if B * B + A * A < 1e-28:
                    continue
                t = 0.25 * np.arctan2(B, A)
                if abs(t) < 1e-13:
                    continue
                c, sn = np.cos(t), np.sin(t)
                for ax in range(4):                     # rotate index `ax`: new_i = c i + s j, new_j = -s i + c j
                    gi = np.take(gb, i, axis=ax).copy()
                    gj = np.take(gb, j, axis=ax).copy()
                    sl_i = [slice(None)] * 4
                    sl_j = [slice(None)] * 4
                    sl_i[ax], sl_j[ax] = i, j
                    gb[tuple(sl_i)] = c * gi + sn * gj
                    gb[tuple(sl_j)] = -sn * gi + c * gj
                ri, rj = R[:, i].copy(), R[:, j].copy()
                R[:, i], R[:, j] = c * ri + sn * rj, -sn * ri + c * rj
                delta = max(delta, abs(t))
        if delta < tol:
            break
    return R


def split_basis_from_integrals(one_body, two_body, nocc: int, order: str = "fiedler"):
    """Integral-only split-localized basis: ER-localize occupied and virtual MOs separately, order by the Fiedler
    vector of the exchange matrix |(pq|qp)| in the localized basis.  Returns (S, occ_mask)."""
    g = np.asarray(two_body, dtype=np.float64)
    n = g.shape[0]
    S = np.zeros((n, n))
    for sl in (np.arange(nocc), np.arange(nocc, n)):
        S[np.ix_(sl, sl)] = _er_localize_block(g, sl)
    occ = np.arange(n) < nocc
    gl = np.einsum("pqrs,pa,qb,rc,sd->abcd", g, S, S, S, S, optimize=True)
    K = np.abs(np.einsum("pqqp->pq", gl))
    np.fill_diagonal(K, 0.0)
    Lap = np.diag(K.sum(1)) - K
    _, ev = np.linalg.eigh(Lap)
    perm = np.argsort(ev[:, 1], kind="stable")
    return S[:, perm], occ[perm]



# ===================================================================================================== public API

class LUCJEnergyTN(LUCJEnergySplitTN):
    """RL reward engine: LUCJ variational energy from a U(1)xU(1) MPS in a split-localized orbital basis.

        ev = LUCJEnergyTN(one_body, two_body, constant, norb, nelec, max_bond=64, device="cuda", name=name)
        E, info = ev.energy(U, Z, t1)          # same state as pretrain.rl.energy.exact_energy(.., make_ucj_op(Z, U,
                                               # "square", t1)); U re-unitarized by its polar factor
    Arguments
      one_body, two_body, constant : active-space integrals (MO basis of the ffsim Hamiltonian, chemists' notation)
      norb, nelec                  : closed shell (nelec[0] == nelec[1])
      max_bond                     : MPS bond dimension chi (accuracy/cost knob, see report)
      device                       : "cuda" (GpuSymMPS, complex64 by default) or "cpu" (NpSymMPS, complex128)
      basis                        : None (default: "boys" if `name` is given, else "er"),
                                     "boys"/"pm" (pyscf localization of occupied and virtual MOs separately, sites
                                     ordered by the Fiedler vector of exp(-centroid distance); needs `name` for
                                     jobs/<name>/<name>.xyz and rhf_dataset/<name>.npz; best measured ordering),
                                     "er" (Edmiston-Ruedenberg from the integrals alone, Fiedler order on exchange
                                     integrals; n29 chi 32: +4.2 mHa vs boys), or an explicit (S, occ_mask) tuple
      zip_margin                   : zip-up bond = zip_margin * max_bond before the optimal canonical compression
    info: discarded weight (sum over all optimal truncations / max single), number of truncations, max bond reached,
          timings (t_state = MPS construction, t_expect = block2 <H>, t_total), imaginary part of <H>.
    """

    def __init__(self, one_body, two_body, constant, norb, nelec, max_bond: int = 64, device: str = "cuda",
                 basis=None, name: str | None = None, dtype=None, cutoff: float | None = None,
                 zip_margin: float = 1.5, mode_tol: float = 1e-8, block2_threads: int = 4,
                 scratch: str | None = None, stack_mem: int = 8 << 30, basis_cache: str | None = None,
                 unitarize: bool = True):
        nocc = int(nelec[0])
        if basis is None:                     # Boys + geometric Fiedler order is ~4 mHa better than ER at n29, chi 32
            basis = "boys" if name is not None else "er"
        if isinstance(basis, tuple):
            S, occ = basis
        elif basis == "er":
            S, occ = split_basis_from_integrals(one_body, two_body, nocc)
        elif basis in ("boys", "pm"):
            assert name is not None, "basis='boys'/'pm' needs the molecule name"
            S, occ = split_localized_basis(name, None, basis, "fiedler", cache_dir=basis_cache)
        else:
            raise ValueError(basis)
        if dtype is None:
            dtype = torch.complex128 if str(device) == "cpu" else torch.complex64
        if cutoff is None:
            cutoff = 1e-12 if dtype == torch.complex128 else 1e-10
        super().__init__(one_body, two_body, constant, norb, nelec, S, occ, max_bond=max_bond, cutoff=cutoff,
                         mode_tol=mode_tol, device=device, dtype=dtype, block2_threads=block2_threads,
                         scratch=scratch, stack_mem=stack_mem, unitarize=unitarize)
        self.zip_margin = zip_margin
