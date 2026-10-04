"""Scalable LUCJ energies with a particle-number-symmetric MPS -- reward engine for RL beyond norb ~17.

Public API
    ev = LUCJEnergyTN(one_body, two_body, constant, norb, nelec, max_bond=128, device="cuda", name=name)
    # basis: "boys" when name is given (geometry), else "er" (integrals only); or an explicit (S, occ_mask)
    E, info = ev.energy(U, Z, t1)          # optional max_bond= per call
    ev.settings()                          # every engine setting -> store it with each result
returns <psi|H|psi> (Hartree, incl. constant) of the ffsim state of pretrain/rl/energy.py
    make_ucj_op(Z, U, "square", t1)  (U re-unitarized by its polar factor, as in the energy jobs)
and info (truncation weights, max bond, timings; see LUCJEnergySplitTN.energy).

Method
  * Orbital basis of the MPS: occupied and virtual MOs localized SEPARATELY (Edmiston-Ruedenberg from the integrals,
    or Boys/PM from the geometry), sites ordered by a Fiedler vector.  |HF> is then an exact product state and the
    LUCJ correlation is local: at norb 15 the final state needs chi(1e-5 discarded per cut) = 159, against 1123 in
    the MO energy order and >= 393 already for |HF> in the network's chemistry-frame chain order.
  * In that basis S the state is  phi_S = exp(i J_1(n^{V_1})) exp(i J_0(n^{V_0})) |HF_S>,  V_k = S^T U_k: every
    square-mask term is a commuting factor exp(i z n_A n_B) = 1 + (e^{iz}-1) b+_A b+_B b_B b_A of two rotated modes,
    an exact bond-17 MPO (Jordan-Wigner strings, U(1)xU(1) labels).  No orbital rotation of the MPS is performed.
  * method="dm" (default on CUDA): every factor is applied with the density-matrix algorithm -- exact Gram
    environments of O|psi> on one side (dense GPU contractions; sector-blocked for bonds >= env_block_min and on the
    CPU), then one sweep that truncates each bond to chi with the exact reduced density matrix rho = M E M^dag: the
    optimal truncation of the exact factor product, no intermediate truncation.  Sweep directions alternate, so no
    canonical moves are needed between factors.  info["discarded_sum"] = sum of the discarded rho weights of all
    truncations, a genuine truncation measure (~ 1 - fidelity to first order).  It is NOT a calibrated energy error,
    but on 32 norb-15..18 tasks (95 runs with err > 0.1 mHa, chi 32-256) err(E) / discarded_sum = 0.9-6.2 Ha (median
    2.2, 10-90 % 1.5-3.7): a usable a-posteriori error bar and a flag for hard cases (raise chi).
  * method="zipup" (legacy; default on the CPU, where the dense environments of "dm" cost far more than the zip-up):
    zip-up with bond <= zip_margin * chi (the Gram of a non-orthonormal right basis -> not an optimal truncation),
    then canonical compression.  Its zip-up truncation weight is reported separately (info["discarded_zip_sum"],
    not Schmidt weights); info["discarded_sum"] then covers the compression sweeps only.  zip_margin is an accuracy
    knob (1.5 -> 3 roughly halves the error at chi 64-128).  With zip_margin 1.5 its error at equal chi was 2-20x
    (median ~5x) that of "dm" on the same 10 tasks (chi 64-256).
  * The t1 rotation and S are absorbed into the Hamiltonian (integrals rotated by final(t1) @ S, real); <H> is the
    block2 quantum-chemistry MPO expectation (SZ, complex MPS, converter verified to 1e-13), on the CPU.
  * Engines: GpuSymMPS (torch; complex64 on CUDA by default, complex128 on CPU) runs both methods; NpSymMPS (numpy
    complex128, CPU) runs "zipup".  Without truncation both reproduce ffsim (complex128: 1e-12 on N2/H2O, ~1e-9 on
    random Hamiltonians with |E| ~ 100 Ha); see pretrain/rl/tests/test_tn_dm.py and test_tn_small.py.

Precision / smoothness: the sector density matrices are diagonalized in fp64 on the host even for complex64 tensors
(fp32 Jacobi on the GPU moved E by 0.7 mHa at chi 64).  complex64 vs complex128 at chi 64/128 (8 tasks, norb 15-17):
|dE| <= 0.09 mHa, mean 0.03.  Identical inputs are bitwise reproducible, but E is only piecewise smooth in (U, Z):
truncation choices switch at near-degenerate weights, so a ~1e-15 change of U (polar factor applied twice) moved the
norb-29 energy by 0.10 / 0.06 / 0.01 mHa at chi 64 / 128 / 256 (norb 16: unchanged).  Do NOT use finite-difference
or autograd gradients through the truncations.

What does not work (kept for the record): TEBD of the Givens-decomposed circuit in the frame chain order
(SymMPS.apply_orbital_rotation / LUCJEnergyTNGivens).  The Clements networks pass through volume-law intermediate
states (norb 15: 0.29 discarded weight for one rotation at chi 128), and the frame-0 -> frame-1 rotation W has
O(1) elements between orbitals 14+ sites apart in the chain at norb 29, so no small-angle/banded path exists.
Fishman-White preparation (fishman_white_sequence) fixes the determinant part but not W.

Validation numbers (accuracy vs exact at norb 15-18, GRPO ranking in 28 perturbation groups, norb-29 chi series,
timings) are produced by pretrain/rl/tests/tn_validate.py + tn_v2_report.py into pretrain/rl/tests/results/tn_v2/.
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
        self.discarded_zip = 0.0       # zip-up Gram truncations (not Schmidt weights), reported separately
        self.max_disc_zip = 0.0
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
                dz = max(0.0, 1.0 - kept_w / total) if total > 0 else 0.0
                self.discarded_zip += dz
                self.max_disc_zip = max(self.max_disc_zip, dz)
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

        ev = LUCJEnergySplitTN(one_body, two_body, constant, norb, nelec, S, occ_mask, max_bond=128, device="cuda")
        E, info = ev.energy(U, Z, t1)
    method "dm" (density-matrix factor application, default) or "zipup" (legacy, zip_margin * chi intermediate).
    """

    def __init__(self, one_body, two_body, constant, norb, nelec, S, occ_mask, max_bond: int = 256,
                 cutoff: float = 1e-12, mode_tol: float = 1e-8, device: str = "cuda", dtype=torch.complex64,
                 block2_threads: int = 8, scratch: str | None = None, stack_mem: int = 8 << 30,
                 unitarize: bool = True, method: str = "dm", zip_margin: float = 2.0):
        import tempfile
        assert nelec[0] == nelec[1], "closed-shell only"
        if method not in ("dm", "zipup"):
            raise ValueError(f"method must be 'dm' or 'zipup', got {method!r}")
        self.h, self.g = np.asarray(one_body, dtype=np.float64), np.asarray(two_body, dtype=np.float64)
        self.const = float(constant)
        self.norb, self.nelec = int(norb), tuple(int(x) for x in nelec)
        self.S = np.asarray(S)
        self.occ = np.asarray(occ_mask, dtype=bool)
        assert self.occ.sum() == self.nelec[0]
        self.max_bond, self.cutoff, self.mode_tol = int(max_bond), float(cutoff), float(mode_tol)
        self.device, self.dtype, self.unitarize = device, dtype, unitarize
        self.method, self.zip_margin = method, float(zip_margin)
        self.basis_desc = "explicit"
        self.block2_threads, self.stack_mem = int(block2_threads), int(stack_mem)
        self._scratch = scratch or tempfile.mkdtemp(prefix="tn_b2_")
        self.b2 = Block2Energy(self.norb, self.nelec, self._scratch, n_threads=block2_threads, stack_mem=stack_mem)
        self._mpo_cache = {}
        self.release_cache = True

    def _engine_cls(self):
        if self.method == "zipup" and str(self.device) == "cpu":
            return NpSymMPS
        return GpuSymMPS

    def settings(self) -> dict:
        """Every setting that affects the energy or its cost (store it with each result)."""
        return {"method": self.method, "zip_margin": self.zip_margin if self.method == "zipup" else None,
                "cutoff": self.cutoff, "dtype": str(self.dtype).replace("torch.", ""), "device": str(self.device),
                "engine": self._engine_cls().__name__, "mode_tol": self.mode_tol, "basis": self.basis_desc,
                "unitarize": self.unitarize, "block2_threads": self.block2_threads,
                "stack_mem_gb": self.stack_mem / (1 << 30), "max_bond_default": self.max_bond}

    def state(self, U, Z, t1=None, max_bond: int | None = None):
        chi = self.max_bond if max_bond is None else int(max_bond)
        U = np.asarray(U, dtype=np.complex128)
        if self.unitarize:
            U = np.stack([polar_unitary(u) for u in U])
        Z = np.asarray(Z, dtype=np.float64)
        occs = [3 if o else 0 for o in self.occ]
        cls = self._engine_cls()
        mps = cls(self.norb, self.nelec, occs) if cls is NpSymMPS else \
            cls(self.norb, self.nelec, occs, device=self.device, dtype=self.dtype)
        for k in range(U.shape[0]):
            V = self.S.conj().T @ U[k]
            facs = []
            for vA, sA, vB, sB, z in square_factors(V, Z[k]):
                facs.append(pair_factor_mpo(vA, sA, vB, sB, z, self.mode_tol))
            facs.sort(key=lambda f: (f[0], f[1]))
            for l, r, Ws, bl, br, dq in facs:
                if self.method == "dm":
                    mps.apply_mpo_dm(l, r, Ws, bl, br, dq, chi, self.cutoff)
                else:
                    mps.apply_mpo_window(l, r, Ws, bl, br, dq, chi, self.cutoff, zip_margin=self.zip_margin)
        if t1 is not None:
            from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
            Fr = orbital_rotation_from_t1_amplitudes(np.asarray(t1, dtype=np.float64)) @ self.S
        else:
            Fr = self.S
        return mps, Fr

    def energy(self, U, Z, t1=None, max_bond: int | None = None):
        """-> (E, info).  info: method, chi, discarded_sum / discarded_max / n_trunc (truncations of the state: for
        "dm" every truncation, a genuine weight; for "zipup" the compression sweeps only), discarded_zip_sum / _max
        (zip-up Gram truncations, "zipup" only, not Schmidt weights), discarded_total, max_bond reached, imag(<H>),
        timings t_state (MPS) / t_env / t_sweep (dm) / t_zip / t_compress (zipup) / t_mpo_build / t_convert /
        t_expect (block2 <H>) / t_total."""
        t0 = time.time()
        chi = self.max_bond if max_bond is None else int(max_bond)
        mps, Fr = self.state(U, Z, t1, chi)
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
        dz = getattr(mps, "discarded_zip", 0.0)
        info = {"method": self.method, "chi": chi, "discarded_sum": mps.discarded, "discarded_max": mps.max_disc,
                "n_trunc": mps.n_trunc, "discarded_zip_sum": dz, "discarded_zip_max": getattr(mps, "max_disc_zip", 0.0),
                "discarded_total": mps.discarded + dz, "max_bond": max(mps.bond_dims()), "imag": e.imag,
                "t_state": t1_ - t0, "t_mpo_build": t2 - t1_, "t_convert": t3 - t2, "t_expect": t4 - t3,
                "t_total": t4 - t0, **{k: v for k, v in mps.stats.items()}}
        del bm, mps
        if str(self.device).startswith("cuda") and self.release_cache:
            torch.cuda.empty_cache()           # several workers share a GPU: return cached blocks after each energy
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
        self.discarded_zip, self.max_disc_zip = 0.0, 0.0
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
                X, newq, kw_, tw_ = self._trunc(M, qrow, qcol, zmax, self.zip_cutoff, account=False)
                dz = max(0.0, 1.0 - kw_ / tw_) if tw_ > 0 else 0.0
                self.discarded_zip += dz
                self.max_disc_zip = max(self.max_disc_zip, dz)
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


def _host_eigh_blocks(Gh, sizes):
    """eigh of the leading sizes[q] x sizes[q] blocks of the stack Gh (numpy), blocks of equal size in one LAPACK batch
    call (padding to size classes was measured slower: it multiplies the O(m^3) work).  Returns [(w ascending, V)]."""
    out = [None] * len(sizes)
    by = {}
    for q, m in enumerate(sizes):
        by.setdefault(int(m), []).append(q)
    for m, qs in by.items():
        S = Gh[qs, :m, :m]
        w, V = np.linalg.eigh(0.5 * (S + S.conj().transpose(0, 2, 1)))
        for i, q in enumerate(qs):
            out[q] = (w[i], V[i])
    return out


class GpuSymMPS(NpSymMPS):
    """SymMPS variant on torch tensors (complex64/128, CUDA or CPU), U(1)xU(1) labels as numpy arrays on the host
    (no device syncs for bookkeeping).  Runs both factor-application methods:
      apply_mpo_dm      density-matrix algorithm (exact environments + one optimal-truncation sweep), default;
      apply_mpo_window  zip-up + canonical compression (legacy).
    Sector blocks of a truncation are gathered into one zero-padded batch (few kernel launches); their Gram /
    density matrices are formed on the device (complex64 or 128) and diagonalized on the host in fp64 (LAPACK,
    equal-size blocks batched), blocks larger than gpu_eig_min with the device eigh.  Fp64 diagonalization of the
    complex64 Grams keeps complex64 within ~0.03 mHa of complex128 (FP32 Jacobi on the device was 0.7 mHa off).
    Bond labels are kept sorted (every truncation / QR returns them in label order)."""

    gpu_eig_min = 256
    batch_eig_max = 0                  # >0: blocks with <= this many rows in one batched device eigh (syevj; FP32
                                       # Jacobi was measured 0.7 mHa off at chi 64 vs fp64 -> off by default)
    env_chunk_elems = 1 << 24          # cap (elements) of the intermediates of one environment step
    env_dense_bytes = 3 << 30          # environments of one window kept dense up to this total size, else label blocks
    env_block_min = 256                # GPU: sector-blocked environment steps once both bonds of a site reach this

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
        self.discarded_zip, self.max_disc_zip = 0.0, 0.0
        self.stats = {"t_zip": 0.0, "t_compress": 0.0, "t_qr": 0.0, "n_mpo": 0, "t_env": 0.0, "t_sweep": 0.0,
                      "n_dm": 0}
        self.zip_cutoff = 1e-15
        self._fast_eig = self.device.type == "cuda" and dtype == torch.complex64
        self.host_eig_dtype = torch.complex128            # host LAPACK in fp64 (also for complex64 Grams)
        # relative floor below which compression singular values are FP noise (zip-up path, split_left)
        self.noise_floor = 1e-7 if dtype == torch.complex64 else 1e-14

    def _t(self, x):
        return torch.as_tensor(x, device=self.device)

    def _sync(self):
        if self.device.type == "cuda":
            torch.cuda.synchronize(self.device)

    # ------------------------------------------------------------------ sector-blocked truncation (batched)
    def _hilbert_cap(self, c, n_left):
        """Largest possible Schmidt rank of label c (left-flow N_alpha * QB + N_beta) at a bond with n_left sites on
        the left: min(dim left sector, dim right sector)."""
        na, nb = divmod(int(c), QB)
        nr = self.n - n_left
        ra, rb = self.nelec[0] - na, self.nelec[1] - nb
        if min(na, nb, ra, rb) < 0 or max(na, nb) > n_left or max(ra, rb) > nr:
            return 0
        return min(math.comb(n_left, na) * math.comb(n_left, nb), math.comb(nr, ra) * math.comb(nr, rb))

    def _padded_blocks(self, M, ridx, cidx):
        """Gather the sector blocks of M into a zero-padded batch: out[t, i, j] = M[ridx[t, i], cidx[t, j]], where
        index len(rows) / len(cols) points at an appended zero row / column.  One transfer, one gather."""
        R_, C_ = M.shape
        Mx = torch.nn.functional.pad(M, (0, 1, 0, 1))
        nb, P = ridx.shape
        K = cidx.shape[1]
        idx = self._t(np.concatenate([ridx.ravel(), cidx.ravel()]))
        return Mx[idx[:nb * P].view(nb, P, 1), idx[nb * P:].view(nb, 1, K)]

    def _env_blocks(self, R, labels, dense=False):
        """Environment matrix R (rows/cols labelled `labels`, the column order of the M it weights) -> metric for
        _trunc: {"dense": R} (kept whole; one gather per truncation), or its label-diagonal blocks in the stable
        label-sorted order {"labels": [c..] (sorted), "pos": {c: t}, "sizes", "E": zero-padded batch (nb, K, K) when
        padding costs < 2x, else None and "blocks": [E_c]} (memory ~ 1/#labels of the dense matrix)."""
        if dense:
            return {"dense": R}
        pc = np.argsort(labels, kind="stable")
        sec = _sectors(labels[pc])
        labs = list(sec)
        sizes = np.array([c1 - c0 for c0, c1 in sec.values()], dtype=np.int64)
        K = int(sizes.max())
        out = {"labels": labs, "pos": {c: t for t, c in enumerate(labs)}, "sizes": sizes, "E": None, "blocks": None}
        if len(labs) * K * K <= 2 * int((sizes ** 2).sum()):
            n = R.shape[0]
            cidx = np.full((len(labs), K), n, dtype=np.int64)
            for t, c in enumerate(labs):
                c0, c1 = sec[c]
                cidx[t, :c1 - c0] = pc[c0:c1]
            out["E"] = self._padded_blocks(R, cidx, cidx)
        else:
            pc_t = self._t(pc)
            out["blocks"] = [R.index_select(0, pc_t[c0:c1]).index_select(1, pc_t[c0:c1]) for c0, c1 in sec.values()]
        return out

    def _trunc(self, M, qrow, qcol, max_bond, cutoff, account=True, metric=None, n_left=None):
        """Dominant row space of the block-diagonal M (rows qrow, cols qcol), optionally with a column metric (an
        _env_blocks dict): rho_c = B_c E_c B_c^dag per label c, else B_c B_c^dag.  All blocks are handled as one
        zero-padded batch (few kernel launches): Grams by batched matmuls; eigen-decompositions of blocks with
        <= batch_eig_max rows by one batched device eigh (complex64 on CUDA), the others on the host (LAPACK).
        Directions beyond the structural rank (#cols of the block; with n_left, the Hilbert-space Schmidt bound of
        the label) are never kept, so cutoff=0 means 'no truncation' without zero-weight bond inflation.
        Returns (X isometry rows x k, labels (k,) sorted, kept weight, total weight)."""
        R_, C_ = M.shape
        pr = np.argsort(qrow, kind="stable")
        pc = np.argsort(qcol, kind="stable")
        secr, secc = _sectors(qrow[pr]), _sectors(qcol[pc])
        labs = [c for c in secr if c in secc]
        nbk = len(labs)
        rsz = np.array([secr[c][1] - secr[c][0] for c in labs], dtype=np.int64)
        csz = np.array([secc[c][1] - secc[c][0] for c in labs], dtype=np.int64)
        P, K = int(rsz.max()), int(csz.max())
        ridx = np.full((nbk, P), R_, dtype=np.int64)
        cidx = np.full((nbk, K), C_, dtype=np.int64)
        for t, c in enumerate(labs):
            r0, r1 = secr[c]
            c0, c1 = secc[c]
            ridx[t, :r1 - r0] = pr[r0:r1]
            cidx[t, :c1 - c0] = pc[c0:c1]
        Bp = self._padded_blocks(M, ridx, cidx)                                  # (nb, P, K)
        if metric is None:
            G = Bp @ Bp.mH
        elif "dense" in metric:
            G = (Bp @ self._padded_blocks(metric["dense"], cidx, cidx)) @ Bp.mH
        elif metric["E"] is not None:
            Ep = metric["E"]
            if len(metric["labels"]) != nbk or metric["labels"] != labs:
                Ep = Ep.index_select(0, self._t(np.array([metric["pos"][c] for c in labs], dtype=np.int64)))
            G = (Bp @ Ep[:, :K, :K]) @ Bp.mH                                     # (nb, P, P)
        else:
            G = torch.zeros((nbk, P, P), dtype=M.dtype, device=self.device)
            for t, c in enumerate(labs):
                m, k = int(rsz[t]), int(csz[t])
                B = Bp[t, :m, :k]
                G[t, :m, :m] = (B @ metric["blocks"][metric["pos"][c]]) @ B.mH
        # eigen-decompositions: small blocks batched on the device (complex64 CUDA), the rest on the host
        small = np.nonzero(rsz <= self.batch_eig_max)[0] if self._fast_eig else np.zeros(0, dtype=np.int64)
        large = np.setdiff1d(np.arange(nbk), small)
        w_blocks = [None] * nbk
        vsrc = [None] * nbk                       # ("dev", t_small) or ("host", V numpy)
        if len(small):
            Ps = int(rsz[small].max())
            Ps = 4 if Ps <= 4 else (8 if Ps <= 8 else (16 if Ps <= 16 else 32))
            Gs = G.index_select(0, self._t(small))[:, :Ps, :Ps]
            Gs = 0.5 * (Gs + Gs.mH)
            pad = np.arange(Ps)[None, :] >= rsz[small][:, None]                # padded rows -> eigenvalue -2 tr
            if pad.any():
                tr = Gs.diagonal(dim1=-2, dim2=-1).real.sum(-1).abs()
                shift = self._t(pad.astype(np.float32)).to(tr.dtype) * (2.0 * tr + 1e-30)[:, None]
                Gs = Gs - torch.diag_embed(shift.to(Gs.dtype))
            ws_dev, Vs_dev = torch.linalg.eigh(Gs)
        else:
            ws_dev, Vs_dev, Ps = None, None, 0
        big = [t for t in large if rsz[t] > self.gpu_eig_min and self.device.type == "cuda"]
        large = np.array([t for t in large if t not in set(big)], dtype=np.int64)
        if len(large):
            Gl = G.index_select(0, self._t(large)).to(self.host_eig_dtype).cpu().numpy()
        if ws_dev is not None:
            wsh = ws_dev.to(torch.float64).cpu().numpy()
            for q, t in enumerate(small):
                m = int(rsz[t])
                w_blocks[t] = wsh[q, Ps - m:]
                vsrc[t] = ("dev", q)
        if len(large):
            for t, (w, V) in zip(large, _host_eigh_blocks(Gl, rsz[large])):
                w_blocks[t] = np.asarray(w, dtype=np.float64)
                vsrc[t] = ("host", V)
        for t in big:                             # very large blocks: device eigh (fallback: host fp64)
            m = int(rsz[t])
            w, V = _herm_eig(G[t, :m, :m])
            w_blocks[t] = w.double().cpu().numpy()
            vsrc[t] = ("devfull", V)
        for t in range(nbk):                      # rare: non-finite device result -> host fp64
            if not np.all(np.isfinite(w_blocks[t])):
                EIG_STATS["fallback"] += 1
                m = int(rsz[t])
                Gq = G[t, :m, :m].detach().to("cpu", torch.complex128).resolve_conj().numpy()
                w, V = np.linalg.eigh(0.5 * (Gq + Gq.conj().T))
                w_blocks[t], vsrc[t] = w, ("host", V)
        ws = []
        for t, c in enumerate(labs):              # rank(B E B^dag) <= #cols: structural zeros are never kept
            w = np.array(w_blocks[t], dtype=np.float64)
            cap = int(csz[t]) if n_left is None else min(int(csz[t]), self._hilbert_cap(c, n_left))
            if len(w) > cap:
                w[:len(w) - cap] = -np.inf
            ws.append(w)
        allw = np.concatenate(ws)
        total = float(np.clip(allw, 0, None).sum())
        order = np.argsort(-allw, kind="stable")
        sw = allw[order]
        n_ok = int(np.isfinite(sw).sum())
        keep = int((sw > cutoff * total).sum()) if cutoff > 0 else n_ok
        keep = max(1, min(keep, max_bond, n_ok))
        kept_w = float(np.clip(sw[:keep], 0, None).sum())
        if account:
            disc = max(0.0, 1.0 - kept_w / total) if total > 0 else 0.0
            self.discarded += disc
            self.max_disc = max(self.max_disc, disc)
            self.n_trunc += 1
        bounds = np.cumsum([0] + [len(w) for w in ws])
        nkeep = np.bincount(np.searchsorted(bounds, order[:keep], side="right") - 1, minlength=nbk)
        # assemble X (rows of M x keep): block t contributes its largest nkeep[t] eigenvectors
        labels = np.empty(keep, dtype=np.int64)
        big_parts = []
        d_src, d_row, d_col = [], [], []
        h_val, h_row, h_col = [], [], []
        j = 0
        for t in range(nbk):
            k = int(nkeep[t])
            if k == 0:
                continue
            m = int(rsz[t])
            rows = ridx[t, :m]
            labels[j:j + k] = labs[t]
            kind, src = vsrc[t]
            if kind == "devfull":
                big_parts.append((rows, j, src[:, m - k:]))
            elif kind == "dev":                   # Vs_dev[q, i, Ps - k + p]  (i < m, p < k)
                q = src
                ii, pp = np.meshgrid(np.arange(m), np.arange(k), indexing="ij")
                d_src.append(((q * Ps + ii) * Ps + (Ps - k + pp)).ravel())
                d_row.append(np.repeat(rows, k))
                d_col.append(np.tile(np.arange(j, j + k), m))
            else:
                h_val.append(src[:, m - k:].ravel())
                h_row.append(np.repeat(rows, k))
                h_col.append(np.tile(np.arange(j, j + k), m))
            j += k
        X = torch.zeros((R_, keep), dtype=M.dtype, device=self.device)
        if d_src:
            ix = self._t(np.concatenate([np.concatenate(d_src), np.concatenate(d_row), np.concatenate(d_col)]))
            nn = ix.numel() // 3
            X.index_put_((ix[nn:2 * nn], ix[2 * nn:]), Vs_dev.reshape(-1)[ix[:nn]].to(M.dtype))
        for rows, j0, Vk in big_parts:
            X[self._t(rows)[:, None], torch.arange(j0, j0 + Vk.shape[1], device=self.device)[None, :]] = Vk.to(M.dtype)
        if h_val:
            hr, hc = np.concatenate(h_row), np.concatenate(h_col)
            ix = self._t(np.concatenate([hr, hc]))
            vals = torch.as_tensor(np.concatenate(h_val), device=self.device).to(M.dtype)
            X.index_put_((ix[:len(hr)], ix[len(hr):]), vals)
        return X, labels, kept_w, total

    # ------------------------------------------------------------------ canonical moves (block QR)
    def _qr_right(self, i):
        """site i -> left isometry (block QR), R absorbed into site i+1."""
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        qrow = (self.q[i][:, None] + SITE_Q_NP[None, :]).reshape(-1)
        qcol = self.q[i + 1]
        pr = np.argsort(qrow, kind="stable")
        pc = np.argsort(qcol, kind="stable")
        secr, secc = _sectors(qrow[pr]), _sectors(qcol[pc])
        pr_t, pc_t = self._t(pr), self._t(pc)
        Ms = A.reshape(Dl * 4, Dr).index_select(0, pr_t).index_select(1, pc_t)
        parts = []
        for c, (c0, c1) in secc.items():
            if c not in secr:
                continue
            r0, r1 = secr[c]
            Q, R = torch.linalg.qr(Ms[r0:r1, c0:c1])
            parts.append((c, r0, r1, c0, c1, Q, R))
        K = sum(p[5].shape[1] for p in parts)
        Qs = torch.zeros((Dl * 4, K), dtype=self.tdtype, device=self.device)
        Rs = torch.zeros((K, Dr), dtype=self.tdtype, device=self.device)
        newq = np.empty(K, dtype=np.int64)
        o = 0
        for c, r0, r1, c0, c1, Q, R in parts:
            k = Q.shape[1]
            Qs[r0:r1, o:o + k] = Q
            Rs[o:o + k, c0:c1] = R
            newq[o:o + k] = c
            o += k
        Qf = torch.empty_like(Qs).index_copy_(0, pr_t, Qs)
        Rf = torch.empty_like(Rs).index_copy_(1, pc_t, Rs)
        self.A[i] = Qf.reshape(Dl, 4, K)
        B = self.A[i + 1]
        self.A[i + 1] = (Rf @ B.reshape(Dr, -1)).reshape(K, 4, B.shape[2])
        self.q[i + 1] = newq
        self.stats["t_qr"] += time.time() - t0

    def _qr_left(self, i):
        """site i -> right isometry (block QR of the conjugate transpose), L absorbed into site i-1."""
        t0 = time.time()
        A = self.A[i]
        Dl, _, Dr = A.shape
        qrow = self.q[i]
        qcol = (self.q[i + 1][None, :] - SITE_Q_NP[:, None]).reshape(-1)
        pr = np.argsort(qrow, kind="stable")
        pc = np.argsort(qcol, kind="stable")
        secr, secc = _sectors(qrow[pr]), _sectors(qcol[pc])
        pr_t, pc_t = self._t(pr), self._t(pc)
        Ms = A.reshape(Dl, 4 * Dr).index_select(0, pr_t).index_select(1, pc_t)
        parts = []
        for c, (r0, r1) in secr.items():
            if c not in secc:
                continue
            c0, c1 = secc[c]
            Q, R = torch.linalg.qr(Ms[r0:r1, c0:c1].mH)          # M_c^dag = Q R  ->  M_c = R^dag Q^dag
            parts.append((c, r0, r1, c0, c1, Q, R))
        K = sum(p[5].shape[1] for p in parts)
        Qs = torch.zeros((K, 4 * Dr), dtype=self.tdtype, device=self.device)      # rows of the new site tensor
        Ls = torch.zeros((Dl, K), dtype=self.tdtype, device=self.device)
        newq = np.empty(K, dtype=np.int64)
        o = 0
        for c, r0, r1, c0, c1, Q, R in parts:
            k = Q.shape[1]
            Qs[o:o + k, c0:c1] = Q.mH
            Ls[r0:r1, o:o + k] = R.mH
            newq[o:o + k] = c
            o += k
        Qf = torch.empty_like(Qs).index_copy_(1, pc_t, Qs)
        Lf = torch.empty_like(Ls).index_copy_(0, pr_t, Ls)
        self.A[i] = Qf.reshape(K, 4, Dr)
        P = self.A[i - 1]
        self.A[i - 1] = (P.reshape(-1, Dl) @ Lf).reshape(P.shape[0], 4, K)
        self.q[i] = newq
        self.stats["t_qr"] += time.time() - t0

    def split_left(self, i, max_bond, cutoff=0.0):
        A = self.A[i]
        Dl, _, Dr = A.shape
        M = A.reshape(Dl, 4 * Dr)
        qcol = (self.q[i + 1][None, :] - SITE_Q_NP[:, None]).reshape(-1)
        Ul, labels, kept_w, total = self._trunc(M, self.q[i], qcol, max_bond, max(cutoff, self.noise_floor))
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
        """[legacy, method="zipup"] zip-up (bond <= zip_margin * max_bond) + canonical compression; ends at l."""
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
                X, newq, kw_, tw_ = self._trunc(M, qrow, qcol, zmax, max(self.zip_cutoff, self.noise_floor),
                                                account=False)
                dz = max(0.0, 1.0 - kw_ / tw_) if tw_ > 0 else 0.0
                self.discarded_zip += dz
                self.max_disc_zip = max(self.max_disc_zip, dz)
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
        self._sync()
        self.stats["t_zip"] += t1 - t0
        self.stats["t_compress"] += time.time() - t1
        self.stats["n_mpo"] += 1

    # ------------------------------------------------------------------ density-matrix application
    def _use_blocked(self, a, b):
        """Sector-blocked environment steps (U(1)xU(1)-sparse A: ~#sectors fewer flops, more launches) on the CPU
        always, on the GPU once the bonds are large enough for the flops to dominate."""
        return self.device.type != "cuda" or min(a, b) >= self.env_block_min

    def _env_step_right(self, A, W, R, qa, qb):
        """R (b*D, b*D) right of site i -> R' (a*D, a*D) right of site i-1; R[x, y] = <Phi(y)|Phi(x)>, x = (b, v).
        qa / qb: (sorted) labels of the bonds left / right of site i (A[a, s, b] = 0 unless qa[a] + q_s = qb[b])."""
        a, _, b = A.shape
        D = W.shape[0]
        Rf = R.reshape(b, D * b * D)
        Wm = W.reshape(D * 4, 4 * D)                                         # [(w,t),(s,v)]
        Wc = W.conj().permute(1, 3, 0, 2).reshape(4 * D, D * 4)             # [(t,V),(x,S)]
        Ac = A.conj().reshape(a, 4 * b)                                      # [A',(S,B)]
        out = torch.empty((a, D, a, D), dtype=self.tdtype, device=self.device)
        step = max(1, int(self.env_chunk_elems // (4 * D * b * D)))
        blocked = self._use_blocked(a, b)
        if blocked:
            seca, secb = _sectors(qa), _sectors(qb)
            pairs = [(s_, seca[beta - SITE_Q_NP[s_]], b0, b1) for beta, (b0, b1) in secb.items() for s_ in range(4)
                     if beta - SITE_Q_NP[s_] in seca]                     # (s, (a0, a1), b0, b1): nonzero A blocks
        for a0 in range(0, a, step):
            ac = min(step, a - a0)
            if blocked:
                X1 = torch.zeros((ac, 4, D * b * D), dtype=self.tdtype, device=self.device)
                for s_, (r0, r1), b0, b1 in pairs:
                    r0c, r1c = max(r0, a0), min(r1, a0 + ac)
                    if r0c < r1c:
                        X1[r0c - a0:r1c - a0, s_] = A[r0c:r1c, s_, b0:b1] @ Rf[b0:b1]
                X1 = X1.reshape(ac, 4 * D, b * D)
            else:
                X1 = (A[a0:a0 + ac].reshape(ac * 4, b) @ Rf).reshape(ac, 4 * D, b * D)    # [a,(s,v),(B,V)]
            X2 = (Wm @ X1).reshape(ac, D, 4, b, D).permute(0, 1, 3, 2, 4).reshape(ac * D * b, 4 * D)
            X3 = (X2 @ Wc).reshape(ac, D, b, D, 4).permute(0, 1, 3, 4, 2).reshape(ac * D * D, 4 * b)
            if blocked:
                X3v = X3.view(ac * D * D, 4, b)
                Rn = torch.zeros((ac * D * D, a), dtype=self.tdtype, device=self.device)
                for s_, (r0, r1), b0, b1 in pairs:
                    Rn[:, r0:r1] += X3v[:, s_, b0:b1] @ Ac.view(a, 4, b)[r0:r1, s_, b0:b1].T
            else:
                Rn = X3 @ Ac.T
            out[a0:a0 + ac] = Rn.reshape(ac, D, D, a).permute(0, 1, 3, 2)                # [a,w,A',x]
        return out.reshape(a * D, a * D)

    def _env_step_left(self, A, W, L, qa, qb):
        """L (a*D, a*D) left of site i -> L' (b*D, b*D) right of site i; L[x, y] = <Phi(y)|Phi(x)>, x = (a, w)."""
        a, _, b = A.shape
        D = W.shape[0]
        Lf = L.reshape(a, D * a * D)                                         # [a,(w,A,W')]
        Wm = W.permute(1, 3, 2, 0).reshape(4 * D, 4 * D)                     # [(t,v),(s,w)]
        Wc = W.conj().permute(1, 0, 3, 2).reshape(4 * D, D * 4)             # [(t,W'),(V,S)]
        Ac = A.conj().reshape(a * 4, b)                                      # [(A,S),B]
        out = torch.empty((b, D, b, D), dtype=self.tdtype, device=self.device)
        step = max(1, int(self.env_chunk_elems // (4 * D * a * D)))
        blocked = self._use_blocked(a, b)
        if blocked:
            seca, secb = _sectors(qa), _sectors(qb)
            pairs = [(s_, r0, r1, secb[alpha + SITE_Q_NP[s_]]) for alpha, (r0, r1) in seca.items() for s_ in range(4)
                     if alpha + SITE_Q_NP[s_] in secb]                    # (s, a0, a1, (b0, b1))
        for c0 in range(0, b, step):
            bc = min(step, b - c0)
            if blocked:
                Y1 = torch.zeros((4, bc, D * a * D), dtype=self.tdtype, device=self.device)
                for s_, r0, r1, (b0, b1) in pairs:
                    b0c, b1c = max(b0, c0), min(b1, c0 + bc)
                    if b0c < b1c:
                        Y1[s_, b0c - c0:b1c - c0] = A[r0:r1, s_, b0c:b1c].T @ Lf[r0:r1]
                Y1 = Y1.reshape(4, bc, D, a * D)
            else:
                Y1 = (A[:, :, c0:c0 + bc].reshape(a, 4 * bc).T @ Lf).reshape(4, bc, D, a * D)   # [s,b,w,(A,W')]
            Y1 = Y1.permute(0, 2, 1, 3).reshape(4 * D, bc * a * D)
            Y2 = (Wm @ Y1).reshape(4, D, bc, a, D).permute(1, 2, 3, 0, 4).reshape(D * bc * a, 4 * D)
            Y3 = (Y2 @ Wc).reshape(D, bc, a, D, 4).permute(0, 1, 3, 2, 4).reshape(D * bc * D, a * 4)
            if blocked:
                Y3v = Y3.view(D * bc * D, a, 4)
                Ln = torch.zeros((D * bc * D, b), dtype=self.tdtype, device=self.device)
                for s_, r0, r1, (b0, b1) in pairs:
                    Ln[:, b0:b1] += Y3v[:, r0:r1, s_] @ Ac.view(a, 4, b)[r0:r1, s_, b0:b1]
            else:
                Ln = Y3 @ Ac
            out[c0:c0 + bc] = Ln.reshape(D, bc, D, b).permute(1, 0, 3, 2)                 # [b,v,B,V]
        return out.reshape(b * D, b * D)

    def _env_dense(self, l, r, D):
        """Keep the window's environments dense if their total size fits env_dense_bytes."""
        item = 8 if self.tdtype == torch.complex64 else 16
        tot = sum((D * len(self.q[i + 1])) ** 2 for i in range(l, r)) * item
        return tot <= self.env_dense_bytes

    def _env_right(self, l, r, Wt, brt, dq):
        """Exact right Gram environments of O|psi> for the bonds right of sites l..r-1, as label blocks."""
        envs = [None] * (r - l)
        if r == l:
            return envs
        D = Wt.shape[1]
        dense = self._env_dense(l, r, D)
        A = self.A[r]
        a, _, b = A.shape
        Wr = torch.einsum("wtsv,v->wts", Wt[r - l], brt)
        Tm = torch.einsum("asb,wts->awtb", A, Wr).reshape(a * D, 4 * b)
        R = Tm @ Tm.mH
        envs[r - 1 - l] = self._env_blocks(R, (self.q[r][:, None] + dq[None, :]).reshape(-1), dense)
        for i in range(r - 1, l, -1):
            R = self._env_step_right(self.A[i], Wt[i - l], R, self.q[i], self.q[i + 1])
            envs[i - 1 - l] = self._env_blocks(R, (self.q[i][:, None] + dq[None, :]).reshape(-1), dense)
        return envs

    def _env_left(self, l, r, Wt, blt, dq):
        """Exact left Gram environments of O|psi> for the bonds right of sites l..r-1, as label blocks."""
        envs = [None] * (r - l)
        if r == l:
            return envs
        D = Wt.shape[1]
        dense = self._env_dense(l, r, D)
        A = self.A[l]
        a, _, b = A.shape
        Wl = torch.einsum("w,wtsv->tsv", blt, Wt[0])
        Tm = torch.einsum("asb,tsv->atbv", A, Wl).reshape(a * 4, b * D)
        L = Tm.T @ Tm.conj()
        envs[0] = self._env_blocks(L, (self.q[l + 1][:, None] + dq[None, :]).reshape(-1), dense)
        for i in range(l + 1, r):
            L = self._env_step_left(self.A[i], Wt[i - l], L, self.q[i], self.q[i + 1])
            envs[i - l] = self._env_blocks(L, (self.q[i + 1][:, None] + dq[None, :]).reshape(-1), dense)
        return envs

    def apply_mpo_dm(self, l, r, Ws, bl, br, dq, max_bond, cutoff=0.0):
        """Density-matrix application of a windowed MPO (identity outside [l, r]): every bond of the window is
        truncated to max_bond (and the relative cutoff) with the exact reduced density matrix of O|psi> (left part
        projected by the previous truncations).  Sites < l must be left- and > r right-isometric: the centre is
        moved into [l, r]; the sweep runs away from the nearer window end (left->right ends at r, else at l)."""
        t0 = time.time()
        c = min(max(self.center, l), r)
        self.move_center(c)
        dq = np.asarray(dq, dtype=np.int64)
        D = len(bl)
        blt = torch.as_tensor(np.asarray(bl), dtype=self.tdtype, device=self.device)
        brt = torch.as_tensor(np.asarray(br), dtype=self.tdtype, device=self.device)
        Wt = torch.as_tensor(np.asarray(Ws), dtype=self.tdtype, device=self.device)
        if (c - l) <= (r - c):                                        # left -> right
            envs = self._env_right(l, r, Wt, brt, dq)
            self._sync()
            t1 = time.time()
            k0 = self.A[l].shape[0]
            C = torch.eye(k0, dtype=self.tdtype, device=self.device)[:, :, None] * blt[None, None, :]   # (k,a,w)
            qleft = self.q[l]
            for i in range(l, r + 1):
                A = self.A[i]
                k, a, _ = C.shape
                b = A.shape[2]
                CA = (C.transpose(1, 2).reshape(k * D, a) @ A.reshape(a, 4 * b)).reshape(k, D, 4, b)
                Tt = torch.tensordot(CA, Wt[i - l], dims=([1, 2], [0, 2]))                  # (k, b, t, v)
                if i < r:
                    M = Tt.permute(0, 2, 1, 3).reshape(k * 4, b * D)
                    qrow = (qleft[:, None] + SITE_Q_NP[None, :]).reshape(-1)
                    qcol = (self.q[i + 1][:, None] + dq[None, :]).reshape(-1)
                    X, newq, _, _ = self._trunc(M, qrow, qcol, max_bond, cutoff, metric=envs[i - l], n_left=i + 1)
                    envs[i - l] = None
                    kn = X.shape[1]
                    self.A[i] = X.reshape(k, 4, kn)
                    C = (X.mH @ M).reshape(kn, b, D)
                    self.q[i + 1] = newq
                    qleft = newq
                else:
                    Tf = torch.tensordot(Tt, brt, dims=([3], [0])).permute(0, 2, 1)
                    self.A[i] = (Tf / torch.linalg.vector_norm(Tf)).contiguous()
            self.center = r
        else:                                                         # right -> left
            envs = self._env_left(l, r, Wt, blt, dq)
            self._sync()
            t1 = time.time()
            k0 = self.A[r].shape[2]
            C = torch.eye(k0, dtype=self.tdtype, device=self.device)[:, None, :] * brt[None, :, None]   # (b,v,k)
            qright = self.q[r + 1]
            for i in range(r, l - 1, -1):
                A = self.A[i]
                W = Wt[i - l]
                a = A.shape[0]
                b, _, k = C.shape
                WC = (W.reshape(D * 16, D) @ C.permute(1, 0, 2).reshape(D, b * k)).reshape(D, 4, 4, b, k)
                WC = WC.permute(2, 3, 0, 1, 4).reshape(4 * b, D * 4 * k)                   # [(s,b),(w,t,k)]
                Tt = (A.reshape(a, 4 * b) @ WC).reshape(a, D, 4, k)                         # [a,w,t,k]
                if i > l:
                    M = Tt.reshape(a * D, 4 * k)
                    qrow = (qright[None, :] - SITE_Q_NP[:, None]).reshape(-1)               # rows of M^T: (t, k)
                    qcol = (self.q[i][:, None] + dq[None, :]).reshape(-1)                   # cols of M^T: (a, w)
                    X, newq, _, _ = self._trunc(M.T, qrow, qcol, max_bond, cutoff, metric=envs[i - 1 - l],
                                                n_left=i)
                    envs[i - 1 - l] = None
                    kn = X.shape[1]
                    self.A[i] = X.T.reshape(kn, 4, k).contiguous()
                    C = (M @ X.conj()).reshape(a, D, kn)
                    self.q[i] = newq
                    qright = newq
                else:
                    Tf = torch.einsum("awtk,w->atk", Tt, blt)
                    self.A[i] = (Tf / torch.linalg.vector_norm(Tf)).contiguous()
            self.center = l
        self._sync()
        t2 = time.time()
        self.stats["t_env"] += t1 - t0
        self.stats["t_sweep"] += t2 - t1
        self.stats["n_dm"] += 1

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

        ev = LUCJEnergyTN(one_body, two_body, constant, norb, nelec, max_bond=128, device="cuda", name=name)
        E, info = ev.energy(U, Z, t1)          # same state as pretrain.rl.energy.exact_energy(.., make_ucj_op(Z, U,
                                               # "square", t1)); U re-unitarized by its polar factor
        ev.settings()                          # every engine setting, to store with the results
    Arguments
      one_body, two_body, constant : active-space integrals (MO basis of the ffsim Hamiltonian, chemists' notation)
      norb, nelec                  : closed shell (nelec[0] == nelec[1])
      max_bond                     : MPS bond dimension chi (accuracy/cost knob; energy(..., max_bond=) overrides)
      method                       : None (default) -> "dm" on CUDA, "zipup" on the CPU.  "dm": density-matrix factor
                                     application, optimal truncation of every exact factor product (GPU: n29 chi 64
                                     ~1 min state build; on one CPU core it was 40x slower than on the GPU).  "zipup":
                                     legacy zip_margin * chi intermediate bond, 2-8x larger error at equal chi.
      device                       : "cuda" (torch engine, complex64 default) or "cpu" (complex128; "dm" runs the
                                     torch engine on the CPU, "zipup" the numpy engine)
      dtype                        : torch.complex64 / torch.complex128 (default: complex64 on CUDA, complex128 on CPU)
      cutoff                       : relative weight below which states are dropped even when chi is not reached
                                     (default 1e-9 for complex64, 1e-12 for complex128); it only removes numerically
                                     empty directions -- chi is the accuracy knob
      basis                        : None (default: "boys" if `name` is given, else "er"),
                                     "boys"/"pm" (pyscf localization of occupied and virtual MOs separately, sites
                                     ordered by the Fiedler vector of exp(-centroid distance); needs `name` for
                                     jobs/<name>/<name>.xyz and rhf_dataset/<name>.npz; best measured ordering),
                                     "er" (Edmiston-Ruedenberg from the integrals alone, Fiedler order on exchange
                                     integrals; n29 chi 32: +4.2 mHa vs boys), or an explicit (S, occ_mask) tuple
      zip_margin                   : method "zipup" only: zip-up bond = zip_margin * max_bond (an accuracy knob:
                                     margin 1.5 -> 3 roughly halves the error at chi 64-128)
      block2_threads               : OpenMP threads of the block2 <H> evaluation (default 1: one core per worker;
                                     t_expect scales ~1/threads)
      stack_mem                    : block2 memory pool in bytes (default 2 GB; virtual reservation per evaluator)
      scratch                      : block2 scratch directory (default: a fresh temp dir); one per worker process
    info: see LUCJEnergySplitTN.energy (truncation weights, max bond, timings).
    """

    def __init__(self, one_body, two_body, constant, norb, nelec, max_bond: int = 128, device: str = "cuda",
                 basis=None, name: str | None = None, dtype=None, cutoff: float | None = None,
                 zip_margin: float = 1.5, mode_tol: float = 1e-8, block2_threads: int = 1,
                 scratch: str | None = None, stack_mem: int = 2 << 30, basis_cache: str | None = None,
                 unitarize: bool = True, method: str | None = None):
        nocc = int(nelec[0])
        if basis is None:                     # Boys + geometric Fiedler order is ~4 mHa better than ER at n29, chi 32
            basis = "boys" if name is not None else "er"
        if isinstance(basis, tuple):
            S, occ = basis
            desc = "explicit"
        elif basis == "er":
            S, occ = split_basis_from_integrals(one_body, two_body, nocc)
            desc = "er/fiedler"
        elif basis in ("boys", "pm"):
            assert name is not None, "basis='boys'/'pm' needs the molecule name"
            S, occ = split_localized_basis(name, None, basis, "fiedler", cache_dir=basis_cache)
            desc = f"{basis}/fiedler"
        else:
            raise ValueError(basis)
        if dtype is None:
            dtype = torch.complex128 if str(device) == "cpu" else torch.complex64
        if cutoff is None:
            cutoff = 1e-12 if dtype == torch.complex128 else 1e-9
        if method is None:                    # dm: dense GPU contractions; on a CPU core its environments cost
            method = "zipup" if str(device) == "cpu" else "dm"     # 10-100x the zip-up (measured), see docstring
        super().__init__(one_body, two_body, constant, norb, nelec, S, occ, max_bond=max_bond, cutoff=cutoff,
                         mode_tol=mode_tol, device=device, dtype=dtype, block2_threads=block2_threads,
                         scratch=scratch, stack_mem=stack_mem, unitarize=unitarize, method=method,
                         zip_margin=zip_margin)
        self.basis_desc = desc
