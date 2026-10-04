"""Exact LUCJ energies on one GPU (torch + numba-CUDA kernels) -- the RL reward for larger molecules.

Reproduces   pretrain.rl.energy.exact_energy(ham, norb, nelec, make_ucj_op(Z, U, conn, t1=t1))
(ffsim statevector <psi|H|psi>) for spin-balanced closed-shell systems (nalpha = nbeta) up to norb 18.

    from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac
    eng = LUCJEnergyGPU.from_npz("C3HN_rxn2586_P")            # rhf_hamiltonians/<name>.npz
    E   = eng.energy(U, Z, t1=t1)                            # U (n_reps,n,n) complex, Z (n_reps,n,n) real
    Es  = eng.energies([(U1, Z1), (U2, Z2)], t1=t1)          # reuses every per-molecule table
    cf  = corr_frac(E, eng.e_hf, eng.e_ccsd)
    python -m pretrain.rl.gpu_energy --tasks <energy_tasks.pkl> --out <json>   # GPU drop-in for the CPU job

State.  psi[a, b] = amplitude of |alpha string a>|beta string b> in ffsim/pyscf order (lexicographic bit
strings; ffsim vectors are vec.reshape(dim_a, dim_b)).  For the spin-balanced UCJ operator on closed-shell
HF, psi is an exactly symmetric matrix (one orbital rotation for both spins, J_aa = J_bb, J_ab symmetric,
HF = e0 e0^T).  Used twice:
  * an orbital rotation psi -> A psi A^T (A = Lambda^k(W) on the string space) is done as a rotation of
    the minor (column) index, an in-place transpose, and the same column rotation:
    (X A^T)^T A^T = A X A^T  because X = X^T;
  * the energy needs only alpha excitation operators (the beta part is the transpose).
  Before the energy psi is symmetrized, (psi + psi^T)/2, which removes the antisymmetric rounding part.

Operator conventions (checked against ffsim 0.0.84, see pretrain/rl/tests/):
  orbital rotation W:  a^dag_i -> sum_j W_ji a^dag_j.  W is factorized on the CPU in float64 into adjacent
    SU(2) Givens blocks and a diagonal D (givens_clustered):  W = M_1 ... M_m D.  On strings, D multiplies
    string S by prod_{p in S} D_p (applied first and merged into the preceding phase pass); M acts on the
    pair (S with p occupied and p+1 empty, S with that electron moved to p+1) as x' = M x and leaves
    strings with both/neither occupied unchanged (det M = 1).  The factorization is ordered so that the
    rotations not touching the top t orbitals come in long runs; a run is applied in ONE pass by a
    shared-memory kernel (strings with a fixed top-t pattern form a contiguous column chunk, and these
    rotations never mix chunks); rotations touching the top orbitals are one pass each.
  UCJ:  W_1 = U_1^dag, W_r = U_r^dag U_{r-1}, W_final = F U_L (F = exp(t1 - t1^dag) via ffsim).
    W_1 acts on HF exactly as a Slater determinant: psi = v v^T, v[a] = det(W_1[occ(a), :k]).
  diagonal Coulomb (time = -1):  psi[a,b] *= exp(i [phi(a) + phi(b) + n_a^T J_ab n_b]),
    phi(s) = sum_{p<q in s} J_aa[p,q] + 1/2 sum_{p in s} J_aa[p,p]  (ffsim uses the upper triangle).
  Hamiltonian (factorize_hamiltonian):  over packed pairs P = (p<=q), F_pq = E_pq + E_qp, F_pp = E_pp,
    H = const + E0 + sum_P g_P Ft_P + 1/2 sum_L s_L (sum_P V_LP Ft_P)^2,   Ft_P = F_P - <HF|F_P|HF>,
    (P|Q) = sum_L s_L V_LP V_LQ (eigendecomposition; |lambda| <= eig_tol * max dropped), h' = h - 1/2 J_rr,
    g = h' + Coulomb(HF), E0 exact in float64.  The HF shift is exact algebra; it shrinks every float32
    quantity by ~25x (complex64 error 1e-5 -> 1e-7 Ha).  With X_L = sum_P V_LP Ft^alpha_P psi:
        E = const + E0 + [2 Re<psi, X_g> + 1/2 sum_L s_L ||X_L + X_L^T||^2] / <psi|psi>.
    X is built tile by tile: rows a of a tile gather their k(n-k)+1 source rows psi[src(a,e), B]
    (index_select) and one batched GEMM with per-row coefficients C_a[L, e] = V_{L,P(a,e)} sign(a,e)
    gives X[A, :, B]; a fused kernel adds the transposed partner tile and reduces (float64 sums).
    Dividing by <psi|psi> (= 1 exactly) removes the first-order effect of complex64 norm drift.

Memory (complex64, dim = C(norb, k)): psi = 8 dim^2 B; coefficient table C = 4 dim (NL+1) (k(n-k)+1) B;
energy tiles 8 T^2 (K + 2(NL+1)) B with T chosen from max_mem_gb; phase-pass blocks ~40 R dim B.
norb 17 (9,8): 4.7 + 1.1 GB (fits a 12 GB card); norb 18 (9,9): 18.9 + 2.5 GB (needs >= ~24 GB).
"""
from __future__ import annotations

import math
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[2]


# ======================================================================= string-space tables (CPU)

@lru_cache(maxsize=None)
def string_tables(norb: int, k: int) -> dict:
    """Tables depending only on (norb, k), cached.  numpy arrays."""
    from pyscf.fci import cistring
    strs = np.asarray(cistring.make_strings(range(norb), k), dtype=np.int64)
    dim = len(strs)
    occ = ((strs[:, None] >> np.arange(norb)[None, :]) & 1).astype(bool)
    occ_list = np.nonzero(occ)[1].reshape(dim, k).astype(np.int64)       # ascending orbitals per string
    # adjacent-orbital pairs: lo = strings with p occ & p+1 empty, hi = same string with p -> p+1
    pairs = []
    for p in range(norb - 1):
        lo = np.nonzero(occ[:, p] & ~occ[:, p + 1])[0]
        hi = cistring.strs2addr(norb, k, strs[lo] ^ ((1 << p) | (1 << (p + 1))))
        pairs.append((lo.astype(np.int32), np.asarray(hi, dtype=np.int32)))
    # packed pairs P = (p <= q), row-major
    pidx = -np.ones((norb, norb), dtype=np.int64)
    P = 0
    for p in range(norb):
        for q in range(p, norb):
            pidx[p, q] = pidx[q, p] = P
            P += 1
    # single excitations seen from the TARGET string a: (F_P psi)[a] = sign * psi[src]
    # pyscf link[a] = (x, y, str1, sign) with E_xy|a> = sign|str1>  =>  <a|E_yx|str1> = sign
    link = cistring.gen_linkstr_index(range(norb), k)
    off = link[:, :, 0] != link[:, :, 1]
    noff = k * (norb - k)
    assert (off.sum(1) == noff).all()
    lo_ = link[off].reshape(dim, noff, 4)
    src = np.concatenate([lo_[:, :, 2], np.arange(dim)[:, None]], axis=1).astype(np.int64)  # last = self
    exc_p = pidx[lo_[:, :, 0], lo_[:, :, 1]].astype(np.int64)
    exc_s = lo_[:, :, 3].astype(np.float64)
    diagP = np.array([pidx[p, p] for p in range(norb)], dtype=np.int64)
    return dict(norb=norb, k=k, dim=dim, strs=strs, occ=occ, occ_list=occ_list, pairs=pairs,
                npair=P, pidx=pidx, src=src, exc_p=exc_p, exc_s=exc_s, diagP=diagP, K=noff + 1)


def givens_clustered(W: np.ndarray, t: int):
    """Exact factorization W = M_1 ... M_m diag(D) (adjacent SU(2) blocks) whose APPLICATION order
    (D, M_m, ..., M_1) consists of segments  LOW_0, TOP_1, LOW_1, ..., TOP_t, LOW_t:
    LOW segments touch only orbitals < n - t (applied by one fused shared-memory pass each),
    TOP_s holds the s rotations that touch the top t orbitals.

    Elimination (A <- G A): columns j = n-1 .. n-t are reduced to e_j by sweeping rows (0,1),
    (1,2), ..., (j-1,j); the remaining (n-t)x(n-t) block is reduced by the plain Givens QR.
    Returns (application list [(p, M)], D)."""
    A = np.array(W, dtype=np.complex128)
    n = A.shape[0]
    t = max(0, min(int(t), n - 1))
    left = []
    for j in range(n - 1, n - t - 1, -1):
        for i in range(j):
            x, y = A[i, j], A[i + 1, j]
            if x == 0:
                continue
            rho = math.hypot(abs(x), abs(y))
            G = np.array([[y, -x], [np.conj(x), np.conj(y)]]) / rho      # G [x, y]^T = [0, rho]
            A[[i, i + 1], :] = G @ A[[i, i + 1], :]
            left.append((i, G))
    m0 = n - t
    for c in range(m0 - 1):
        for r in range(m0 - 1, c, -1):
            x, y = A[r - 1, c], A[r, c]
            if y == 0:
                continue
            rho = math.hypot(abs(x), abs(y))
            G = np.array([[np.conj(x), np.conj(y)], [-y, x]]) / rho      # G [x, y]^T = [rho, 0]
            A[[r - 1, r], :] = G @ A[[r - 1, r], :]
            left.append((r - 1, G))
    return [(p, G.conj().T) for p, G in reversed(left)], np.diag(A).copy()


def split_segments(apply_list, n: int, t: int):
    """[(kind, [(p, M), ...])], kind 'low' if p + 1 < n - t else 'top'; consecutive runs merged."""
    segs = []
    for p, M in apply_list:
        kind = "low" if p + 1 < n - t else "top"
        if segs and segs[-1][0] == kind:
            segs[-1][1].append((p, M))
        else:
            segs.append((kind, [(p, M)]))
    return segs


@lru_cache(maxsize=None)
def _pairs_only(norb: int, k: int):
    """Adjacent-orbital pair tables of the (norb, k) string space (k may be 0 or norb)."""
    from pyscf.fci import cistring
    if k < 0 or k > norb:
        return np.zeros(0, np.int64), []
    strs = np.asarray(cistring.make_strings(range(norb), k), dtype=np.int64) if norb > 0 else np.zeros(1, np.int64)
    occ = ((strs[:, None] >> np.arange(norb)[None, :]) & 1).astype(bool)
    pairs = []
    for p in range(norb - 1):
        lo = np.nonzero(occ[:, p] & ~occ[:, p + 1])[0]
        hi = np.searchsorted(strs, strs[lo] ^ ((1 << p) | (1 << (p + 1))))
        pairs.append((lo.astype(np.int32), hi.astype(np.int32)))
    return strs, pairs


@lru_cache(maxsize=None)
def chunk_tables(norb: int, k: int, t: int) -> dict:
    """Row chunks for the fused low-segment kernel: strings with equal occupation of the top t orbitals
    form a contiguous row range (lexicographic order) whose low part is the canonical (norb-t, k-h)
    string space.  Local adjacent-pair tables per k_low = k - h."""
    tb = string_tables(norb, k)
    strs = tb["strs"]
    nl = norb - t
    top = strs >> nl
    starts = np.r_[0, np.nonzero(np.diff(top))[0] + 1]
    lens = np.diff(np.r_[starts, len(strs)])
    kls = np.array([k - bin(int(top[s])).count("1") for s in starts], dtype=np.int32)
    kl_vals = sorted(set(kls.tolist()))
    npl = max(nl - 1, 1)
    off = np.zeros((k + 1, npl), dtype=np.int32)
    cnt = np.zeros((k + 1, npl), dtype=np.int32)
    lo_all, hi_all, pos = [], [], 0
    for kl in kl_vals:
        _, pairs = _pairs_only(nl, kl)
        for p, (lo, hi) in enumerate(pairs):
            off[kl, p] = pos
            cnt[kl, p] = len(lo)
            lo_all.append(lo)
            hi_all.append(hi)
            pos += len(lo)
    cat = (lambda xs: np.concatenate(xs).astype(np.int32) if xs else np.zeros(1, np.int32))
    lo, hi = cat(lo_all), cat(hi_all)
    # per-k_low blocks of the flat pair arrays (staged into shared memory by the minor-mode kernel)
    kl_base = np.zeros(k + 1, dtype=np.int32)
    kl_len = np.zeros(k + 1, dtype=np.int32)
    for kl in kl_vals:
        kl_base[kl] = off[kl, 0]
        kl_len[kl] = int(cnt[kl].sum())
    off_rel = (off - kl_base[:, None]).astype(np.int32)
    assert int(lens.max()) < 65536
    return dict(t=t, starts=starts.astype(np.int32), lens=lens.astype(np.int32), kls=kls, off=off, cnt=cnt,
                lo=lo, hi=hi, lo16=lo.astype(np.uint16), hi16=hi.astype(np.uint16), kl_base=kl_base,
                kl_len=kl_len, off_rel=off_rel, max_len=int(lens.max()), max_kl_len=int(kl_len.max()))


def plan_fusion(norb: int, k: int, itemsize: int, smem_bytes: int = 48 * 1024, min_rows: int = 2,
                max_rows: int = 16):
    """(t, R) for the minor-mode fused kernel: smallest t such that R >= min_rows rows of the largest chunk
    plus the chunk's uint16 pair tables fit in smem_bytes.  None if no t works."""
    for t in range(1, norb - 1):
        ct = chunk_tables(norb, k, t)
        if ct["max_len"] >= 65536:
            continue
        R = (smem_bytes - 4 * ct["max_kl_len"] - 16) // (ct["max_len"] * itemsize)
        if R >= min_rows:
            return t, int(min(R, max_rows))
    return None


def factorize_hamiltonian(one_body, two_body, k: int, eig_tol: float = 1e-12):
    """HF-shifted factorized Hamiltonian.

    H = const + E0 + sum_P g_P Ft_P + 1/2 sum_L s_L (sum_P V_LP Ft_P)^2,   Ft_P = F_P - <HF|F_P|HF>,
    with (P|Q) = sum_L s_L V_LP V_LQ (s_L = +1 for L < npos, -1 after), g = h' + sum_L s_L v_L V_L,
    v_L = sum_P V_LP n_P, n_P = <HF|F_P|HF> (= 2 for P = (p,p), p < k), E0 = h'.n + 1/2 sum_L s_L v_L^2.
    The shift is exact algebra; it makes every float32 quantity correlation-sized (the HF determinant
    contributes exactly zero to the GEMM part).  Returns (g (P,), V (NL, P), npos, E0, info)."""
    h = np.asarray(one_body, dtype=np.float64)
    eri = np.asarray(two_body, dtype=np.float64)
    n = h.shape[0]
    hp = h - 0.5 * np.einsum("prrq->pq", eri)
    iu = np.triu_indices(n)
    Wm = eri[iu[0], iu[1]][:, iu[0], iu[1]]          # (P|Q), P, Q row-major packed p<=q
    asym = float(np.abs(Wm - Wm.T).max())
    Wm = 0.5 * (Wm + Wm.T)
    lam, vec = np.linalg.eigh(Wm)
    thr = eig_tol * max(np.abs(lam).max(), 1e-300)
    pos = lam > thr
    neg = lam < -thr
    V = np.concatenate([(vec[:, pos] * np.sqrt(lam[pos])).T, (vec[:, neg] * np.sqrt(-lam[neg])).T], 0)
    s = np.concatenate([np.ones(int(pos.sum())), -np.ones(int(neg.sum()))])
    nP = np.zeros(len(iu[0]))
    nP[(iu[0] == iu[1]) & (iu[0] < k)] = 2.0
    hpP = hp[iu]
    v = V @ nP
    g = hpP + (s * v) @ V
    E0 = float(hpP @ nP + 0.5 * np.sum(s * v * v))
    dropped = lam[~(pos | neg)]
    info = dict(n_pairs=len(lam), n_pos=int(pos.sum()), n_neg=int(neg.sum()),
                lam_max=float(lam.max()), lam_min=float(lam.min()),
                dropped_abs_sum=float(np.abs(dropped).sum()), eri_asym=asym, E0=E0)
    return g, V, int(pos.sum()), E0, info


def corr_frac(E, e_hf, e_ccsd):
    """Fraction of the CCSD correlation energy recovered: (E_HF - E) / (E_HF - E_CCSD)."""
    return (e_hf - E) / (e_hf - e_ccsd)


def interaction_masks(connectivity: str, norb: int):
    from pretrain.rl.energy import interaction_pairs
    out = []
    for pairs in interaction_pairs(connectivity, norb):
        if pairs is None:
            out.append(np.ones((norb, norb), dtype=bool))
            continue
        m = np.zeros((norb, norb), dtype=bool)
        for i, j in pairs:
            m[i, j] = m[j, i] = True
        out.append(m)
    return out


# ======================================================================= GPU kernels

_NB_KERNELS: dict = {}


def _numba_kernels(cplx: str):
    """numba-CUDA kernels for complex64 ('c8') or complex128 ('c16') states."""
    if cplx in _NB_KERNELS:
        return _NB_KERNELS[cplx]
    from numba import cuda, types

    nbc = types.complex64 if cplx == "c8" else types.complex128
    nbr = types.float32 if cplx == "c8" else types.float64
    f64 = types.float64
    u32 = types.uint32
    u16 = types.uint16

    @cuda.jit
    def givens_rows(X, lo, hi, m00, m01, m10, m11):
        c = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
        p = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
        if p < lo.shape[0] and c < X.shape[1]:
            i = lo[p]
            j = hi[p]
            x = X[i, c]
            y = X[j, c]
            X[i, c] = m00 * x + m01 * y
            X[j, c] = m10 * x + m11 * y

    @cuda.jit
    def givens_cols(X, lo, hi, m00, m01, m10, m11):
        # one adjacent-orbital rotation on the MINOR (column) index of every row
        q = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
        r = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
        if q < lo.shape[0] and r < X.shape[0]:
            i = lo[q]
            j = hi[q]
            x = X[r, i]
            y = X[r, j]
            X[r, i] = m00 * x + m01 * y
            X[r, j] = m10 * x + m11 * y

    @cuda.jit
    def givens_low_minor(X, ch_start, ch_len, ch_kl, rot_p, rot_m, plo, phi, kl_base, kl_len, off_rel, cnt_tab,
                         R, stride, u_off):
        # One block = R consecutive rows x one chunk of the MINOR index (strings with a fixed top-orbital
        # pattern = a contiguous column range), held in shared memory while the whole low segment of
        # Givens rotations is applied.  The chunk's local pair tables (uint16) are staged in shared memory
        # behind the data.  Loads/stores are contiguous row segments.  uint32 index arithmetic.
        Sc = cuda.shared.array(0, nbc)
        Su = cuda.shared.array(0, u16)
        c = cuda.blockIdx.y
        Ru = u32(R)
        st = u32(stride)
        b0 = u32(cuda.blockIdx.x) * Ru
        c0 = u32(ch_start[c])
        L = u32(ch_len[c])
        kl = ch_kl[c]
        nrow = u32(X.shape[0])
        tid = u32(cuda.threadIdx.x)
        nt = u32(cuda.blockDim.x)
        kb = u32(kl_base[kl])
        kn = u32(kl_len[kl])
        uo = u32(u_off)
        i = tid
        while i < kn:
            Su[uo + i] = plo[kb + i]
            Su[uo + kn + i] = phi[kb + i]
            i += nt
        for r in range(R):
            b = b0 + u32(r)
            if b < nrow:
                base = u32(r) * st
                i = tid
                while i < L:
                    Sc[base + i] = X[b, c0 + i]
                    i += nt
        cuda.syncthreads()
        for rr in range(rot_p.shape[0]):
            p = rot_p[rr]
            m00 = rot_m[rr, 0]
            m01 = rot_m[rr, 1]
            m10 = rot_m[rr, 2]
            m11 = rot_m[rr, 3]
            o = uo + u32(off_rel[kl, p])
            cnt = u32(cnt_tab[kl, p])
            q = tid
            while q < cnt:
                i1 = u32(Su[o + q])
                i2 = u32(Su[o + kn + q])
                for r in range(R):
                    base = u32(r) * st
                    x = Sc[base + i1]
                    y = Sc[base + i2]
                    Sc[base + i1] = m00 * x + m01 * y
                    Sc[base + i2] = m10 * x + m11 * y
                q += nt
            cuda.syncthreads()
        for r in range(R):
            b = b0 + u32(r)
            if b < nrow:
                base = u32(r) * st
                i = tid
                while i < L:
                    X[b, c0 + i] = Sc[base + i]
                    i += nt

    @cuda.jit
    def phase_diag(X, strs, rowfac, ez, n, init):
        # diagonal-Coulomb phase for diagonal J_ab (square/hex/heavy-hex), complex128 arithmetic:
        # X[r, c] = (init ? 1 : X[r, c]) * rowfac[r] rowfac[c] prod_{p in s_r & s_c} exp(i J_ab[p, p])
        # rows on grid x (gridDim.y is capped at 65535 < dim for norb 19)
        c = cuda.blockIdx.y * cuda.blockDim.x + cuda.threadIdx.x
        r = cuda.blockIdx.x
        if c < X.shape[1]:
            m = strs[r] & strs[c]
            ph = rowfac[r] * rowfac[c]
            for p in range(n):
                if (m >> p) & 1:
                    ph *= ez[p]
            if init:
                X[r, c] = ph
            else:
                X[r, c] = X[r, c] * ph

    @cuda.jit
    def transpose_sq(X):
        t1 = cuda.shared.array((32, 33), nbc)
        t2 = cuda.shared.array((32, 33), nbc)
        bj = cuda.blockIdx.x
        bi = cuda.blockIdx.y
        if bj < bi:
            return
        n = X.shape[0]
        tx = cuda.threadIdx.x
        ty = cuda.threadIdx.y
        for k in range(0, 32, 8):
            r = bi * 32 + ty + k
            c = bj * 32 + tx
            if r < n and c < n:
                t1[ty + k, tx] = X[r, c]
            r = bj * 32 + ty + k
            c = bi * 32 + tx
            if r < n and c < n:
                t2[ty + k, tx] = X[r, c]
        cuda.syncthreads()
        for k in range(0, 32, 8):
            r = bi * 32 + ty + k
            c = bj * 32 + tx
            if r < n and c < n:
                X[r, c] = t2[tx, ty + k]
            r = bj * 32 + ty + k
            c = bi * 32 + tx
            if r < n and c < n:
                X[r, c] = t1[tx, ty + k]

    @cuda.jit
    def symmetrize_sq(X):
        # X <- (X + X^T) / 2 in place (projects out the antisymmetric rounding part)
        t1 = cuda.shared.array((32, 33), nbc)
        t2 = cuda.shared.array((32, 33), nbc)
        bj = cuda.blockIdx.x
        bi = cuda.blockIdx.y
        if bj < bi:
            return
        n = X.shape[0]
        tx = cuda.threadIdx.x
        ty = cuda.threadIdx.y
        half = nbr(0.5)
        for k in range(0, 32, 8):
            r = bi * 32 + ty + k
            c = bj * 32 + tx
            if r < n and c < n:
                t1[ty + k, tx] = X[r, c]
            r = bj * 32 + ty + k
            c = bi * 32 + tx
            if r < n and c < n:
                t2[ty + k, tx] = X[r, c]
        cuda.syncthreads()
        for k in range(0, 32, 8):
            r = bi * 32 + ty + k
            c = bj * 32 + tx
            if r < n and c < n:
                X[r, c] = (t1[ty + k, tx] + t2[tx, ty + k]) * half
            r = bj * 32 + ty + k
            c = bi * 32 + tx
            if r < n and c < n:
                X[r, c] = (t2[ty + k, tx] + t1[tx, ty + k]) * half

    @cuda.jit
    def epilogue(X1, X2, P1, P2, nq, npos, diag, out):
        # X1 = X[A, :, B] (nA, NL+1, nB), X2 = X[B, :, A]; P1 = psi[A, B], P2 = psi[B, A]
        # out[block] = (sum_{L<npos} |Y_L|^2, sum_{npos<=L<nq} |Y_L|^2, Re<P1,X1_h> (+ Re<P2,X2_h>), |P|^2)
        # with Y_L[a, b] = X1[a, L, b] + X2[b, L, a] over the (A, B) tile.
        tile = cuda.shared.array((32, 33), nbc)
        red = cuda.shared.array((4, 256), f64)
        tx = cuda.threadIdx.x
        ty = cuda.threadIdx.y
        a0 = cuda.blockIdx.y * 32
        b0 = cuda.blockIdx.x * 32
        nA = X1.shape[0]
        nB = X1.shape[2]
        h = X1.shape[1] - 1
        s_pos = 0.0
        s_neg = 0.0
        for L in range(nq):
            for k in range(0, 32, 8):
                b = b0 + ty + k
                a = a0 + tx
                if b < nB and a < nA:
                    tile[ty + k, tx] = X2[b, L, a]
            cuda.syncthreads()
            part = nbr(0.0)
            for k in range(0, 32, 8):
                a = a0 + ty + k
                b = b0 + tx
                if a < nA and b < nB:
                    y = X1[a, L, b] + tile[tx, ty + k]
                    part += y.real * y.real + y.imag * y.imag
            if L < npos:
                s_pos += part
            else:
                s_neg += part
            cuda.syncthreads()
        e1 = 0.0
        nr = 0.0
        for k in range(0, 32, 8):
            a = a0 + ty + k
            b = b0 + tx
            if a < nA and b < nB:
                p = P1[a, b]
                x = X1[a, h, b]
                e1 += f64(p.real) * f64(x.real) + f64(p.imag) * f64(x.imag)
                nr += f64(p.real) * f64(p.real) + f64(p.imag) * f64(p.imag)
            if diag == 0:
                b = b0 + ty + k
                a = a0 + tx
                if b < nB and a < nA:
                    p = P2[b, a]
                    x = X2[b, h, a]
                    e1 += f64(p.real) * f64(x.real) + f64(p.imag) * f64(x.imag)
                    nr += f64(p.real) * f64(p.real) + f64(p.imag) * f64(p.imag)
        t = ty * 32 + tx
        red[0, t] = s_pos
        red[1, t] = s_neg
        red[2, t] = e1
        red[3, t] = nr
        cuda.syncthreads()
        s = 128
        while s > 0:
            if t < s:
                for q in range(4):
                    red[q, t] += red[q, t + s]
            cuda.syncthreads()
            s //= 2
        if t == 0:
            bid = cuda.blockIdx.y * cuda.gridDim.x + cuda.blockIdx.x
            for q in range(4):
                out[bid, q] = red[q, 0]

    _NB_KERNELS[cplx] = dict(givens_rows=givens_rows, givens_cols=givens_cols, givens_low_minor=givens_low_minor,
                             phase_diag=phase_diag,
                             transpose_sq=transpose_sq, symmetrize_sq=symmetrize_sq, epilogue=epilogue)
    return _NB_KERNELS[cplx]


class _NumbaBackend:
    name = "numba"

    def __init__(self, device: torch.device, cdtype: torch.dtype):
        from numba import cuda
        self.cuda = cuda
        self.cplx = "c8" if cdtype == torch.complex64 else "c16"
        self.k = _numba_kernels(self.cplx)
        self.np_c = np.complex64 if cdtype == torch.complex64 else np.complex128
        self.dev_index = device.index if device.index is not None else torch.cuda.current_device()

    def _stream(self):
        return self.cuda.external_stream(torch.cuda.current_stream().cuda_stream)

    def arr(self, t: torch.Tensor):
        return self.cuda.as_cuda_array(t)

    def givens_sequence(self, psi_nb, ncol: int, rots, pair_nb):
        """Apply [(p, M)] in the given order to the rows of psi (one pass per rotation)."""
        st = self._stream()
        kern = self.k["givens_rows"]
        c = self.np_c
        with self.cuda.gpus[self.dev_index]:
            for p, M in rots:
                lo, hi = pair_nb[p]
                m = lo.shape[0]
                grid = ((ncol + 127) // 128, (m + 1) // 2)
                kern[grid, (128, 2), st](psi_nb, lo, hi, c(M[0, 0]), c(M[0, 1]), c(M[1, 0]), c(M[1, 1]))

    def upload_segments(self, segments):
        """One host->device copy of every fused segment's rotation table: [(kind, rots, dev_p, dev_m)].
        (A pageable copy per segment would synchronize the stream before each fused pass.)"""
        lows = [rots for kind, rots in segments if kind == "low"]
        if not lows:
            return [(kind, rots, None, None) for kind, rots in segments]
        P = np.array([p for rots in lows for p, _ in rots], dtype=np.int32)
        Mh = np.array([[M[0, 0], M[0, 1], M[1, 0], M[1, 1]] for rots in lows for _, M in rots], dtype=self.np_c)
        with self.cuda.gpus[self.dev_index]:
            dP = self.cuda.to_device(P, stream=self._stream())
            dM = self.cuda.to_device(Mh, stream=self._stream())
        out, o = [], 0
        for kind, rots in segments:
            if kind == "low":
                out.append((kind, rots, dP[o:o + len(rots)], dM[o:o + len(rots)]))
                o += len(rots)
            else:
                out.append((kind, rots, None, None))
        return out

    def givens_minor(self, psi_nb, nrow: int, segments, pair_nb, ch):
        """Apply [(kind, rots[, dev_p, dev_m])] to the MINOR (column) index of every row: 'low' segments in
        one fused shared-memory pass over contiguous column chunks, 'top' rotations one pass each."""
        st = self._stream()
        kl_ = self.k["givens_low_minor"]
        kc = self.k["givens_cols"]
        c = self.np_c
        if segments and len(segments[0]) == 2:
            segments = self.upload_segments(segments)
        with self.cuda.gpus[self.dev_index]:
            for kind, rots, rp, rm in segments:
                if kind == "top":
                    for p, M in rots:
                        lo, hi = pair_nb[p]
                        grid = ((lo.shape[0] + 127) // 128, (nrow + 1) // 2)
                        kc[grid, (128, 2), st](psi_nb, lo, hi, c(M[0, 0]), c(M[0, 1]), c(M[1, 0]), c(M[1, 1]))
                    continue
                grid = ((nrow + ch["R"] - 1) // ch["R"], ch["nch"])
                kl_[grid, 256, st, ch["smem_minor"]](psi_nb, ch["starts"], ch["lens"], ch["kls"], rp, rm,
                                                     ch["lo16"], ch["hi16"], ch["kl_base"], ch["kl_len"],
                                                     ch["off_rel"], ch["cnt"], ch["R"], ch["stride"], ch["u_off"])

    def phase_diag(self, psi_nb, strs_nb, rowfac, ez, n: int, init: bool):
        st = self._stream()
        dim = psi_nb.shape[1]
        with self.cuda.gpus[self.dev_index]:
            self.k["phase_diag"][(psi_nb.shape[0], (dim + 255) // 256), 256, st](
                psi_nb, strs_nb, self.arr(rowfac), self.arr(ez), n, 1 if init else 0)

    def transpose(self, psi_nb, n: int):
        st = self._stream()
        nt = (n + 31) // 32
        with self.cuda.gpus[self.dev_index]:
            self.k["transpose_sq"][(nt, nt), (32, 8), st](psi_nb)

    def symmetrize(self, psi_nb, n: int):
        st = self._stream()
        nt = (n + 31) // 32
        with self.cuda.gpus[self.dev_index]:
            self.k["symmetrize_sq"][(nt, nt), (32, 8), st](psi_nb)

    def epilogue(self, X1, X2, P1, P2, nq, npos, diag, out):
        st = self._stream()
        nA, nB = X1.shape[0], X1.shape[2]
        grid = ((nB + 31) // 32, (nA + 31) // 32)
        with self.cuda.gpus[self.dev_index]:
            self.k["epilogue"][grid, (32, 8), st](self.arr(X1), self.arr(X2), self.arr(P1), self.arr(P2),
                                                  nq, npos, 1 if diag else 0, self.arr(out))
        return grid[0] * grid[1]


class _TorchBackend:
    """Pure-torch fallback (no numba): same math, more memory traffic."""
    name = "torch"

    def __init__(self, device, cdtype, col_chunk: int = 4096):
        self.col_chunk = col_chunk

    def arr(self, t):
        return t

    def givens_sequence(self, psi, ncol, rots, pairs):
        for p, M in rots:
            lo, hi = pairs[p]
            for c0 in range(0, ncol, self.col_chunk):
                v = psi[:, c0:c0 + self.col_chunk]
                x = v.index_select(0, lo)
                y = v.index_select(0, hi)
                v.index_copy_(0, lo, x * complex(M[0, 0]) + y * complex(M[0, 1]))
                v.index_copy_(0, hi, x * complex(M[1, 0]) + y * complex(M[1, 1]))

    def transpose(self, psi, n, B: int = 2048):
        for i in range(0, n, B):
            for j in range(i, n, B):
                if i == j:
                    psi[i:i + B, i:i + B] = psi[i:i + B, i:i + B].t().clone()
                else:
                    tmp = psi[i:i + B, j:j + B].clone()
                    psi[i:i + B, j:j + B] = psi[j:j + B, i:i + B].t()
                    psi[j:j + B, i:i + B] = tmp.t()

    def symmetrize(self, psi, n, B: int = 2048):
        for i in range(0, n, B):
            for j in range(i, n, B):
                s = 0.5 * (psi[i:i + B, j:j + B] + psi[j:j + B, i:i + B].t())
                psi[i:i + B, j:j + B] = s
                psi[j:j + B, i:i + B] = s.t()

    def epilogue(self, X1, X2, P1, P2, nq, npos, diag, out, chunk: int = 8):
        # chunks over L keep the temporaries small; float64 accumulation
        s_pos = torch.zeros((), device=X1.device, dtype=torch.float64)
        s_neg = torch.zeros((), device=X1.device, dtype=torch.float64)
        for l0 in range(0, nq, chunk):
            l1 = min(nq, l0 + chunk)
            Y = torch.view_as_real(X1[:, l0:l1] + X2[:, l0:l1].permute(2, 1, 0))
            if l1 <= npos:
                s_pos += torch.sum(Y * Y, dtype=torch.float64)
            elif l0 >= npos:
                s_neg += torch.sum(Y * Y, dtype=torch.float64)
            else:
                sq = torch.sum(Y * Y, dim=(0, 2, 3), dtype=torch.float64)
                s_pos += sq[: npos - l0].sum()
                s_neg += sq[npos - l0:].sum()
            del Y
        pr1, xr1 = torch.view_as_real(P1), torch.view_as_real(X1[:, -1])
        e1 = torch.sum(pr1 * xr1, dtype=torch.float64)
        nr = torch.sum(pr1 * pr1, dtype=torch.float64)
        if not diag:
            pr2, xr2 = torch.view_as_real(P2), torch.view_as_real(X2[:, -1])
            e1 = e1 + torch.sum(pr2 * xr2, dtype=torch.float64)
            nr = nr + torch.sum(pr2 * pr2, dtype=torch.float64)
        out[0, 0] = s_pos
        out[0, 1] = s_neg
        out[0, 2] = e1
        out[0, 3] = nr
        return 1


# ======================================================================= engine

def _default_ham_dir() -> Path:
    return ROOT / "rhf_hamiltonians"


class LUCJEnergyGPU:
    """Exact <psi|H|psi> of the spin-balanced LUCJ state on one GPU.

    Args:
      one_body (n,n), two_body (n,n,n,n) real (chemist notation, ffsim MolecularHamiltonian), constant.
      norb, nelec = (k, k).
      dtype: torch.complex64 (default; |dE| <~ 1e-7 Ha vs ffsim on the stored LUCJ tasks) or
        torch.complex128 (reference; fp64 GEMMs, slow on consumer GPUs).
      max_mem_gb: GPU memory the engine may use (state + tables + workspace); default: free memory - 0.6 GB
        (the energy tiles take what is left; larger tiles are a little faster).
      eig_tol: drop ERI-matrix eigenvalues below eig_tol * max.  1e-12 (default) is exact to ~1e-12;
        1e-10 changes E by <~1e-9 Ha; 1e-8 removes ~15-20% of the GEMM rows at a ~1e-7 Ha downward bias.
      polar: re-unitarize U by its SVD polar factor before use (what the CPU energy job does).
      backend: "numba" (custom kernels), "torch" (fallback, ~4-5x slower) or "auto".
      fuse / smem_kb / fuse_t / fuse_rows: fused shared-memory Givens passes (numba) and their sizing.
      max_tile: cap of the energy tile (1024; larger tiles were slower on the TITAN Xp; retune per GPU).
      alloc_state=False, tile: tables only / fixed tile, for component benchmarks on partial states.
    Methods: energy, energies (batched over parameter sets), state, energy_of_state, set_hamiltonian /
      set_npz (switch molecule, same (norb, nelec), keeps the state buffer), release, corr_frac.
    """

    def __init__(self, one_body, two_body, constant, norb: int, nelec, device="cuda",
                 dtype=torch.complex64, max_mem_gb: float | None = None, eig_tol: float = 1e-12,
                 polar: bool = True, backend: str = "auto", e_hf: float | None = None,
                 e_ccsd: float | None = None, name: str | None = None, fuse: bool = True,
                 smem_kb: float = 48, fuse_t: int | None = None, fuse_rows: int | None = None,
                 alloc_state: bool = True, tile: int | None = None, max_tile: int = 1024):
        nelec = tuple(int(x) for x in nelec)
        if nelec[0] != nelec[1]:
            raise ValueError(f"closed-shell (nalpha == nbeta) only, got {nelec}")
        self.norb, self.k, self.nelec = int(norb), nelec[0], nelec
        self.device = torch.device(device)
        if self.device.type != "cuda":
            raise ValueError("LUCJEnergyGPU needs a CUDA device")
        if self.device.index is None:
            self.device = torch.device("cuda", torch.cuda.current_device())
        if dtype not in (torch.complex64, torch.complex128):
            raise ValueError(dtype)
        self.cdtype = dtype
        self.rdtype = torch.float32 if dtype == torch.complex64 else torch.float64
        self.constant = float(constant)
        self.polar = polar
        self.e_hf, self.e_ccsd, self.name = e_hf, e_ccsd, name
        self.timing: dict = {}

        tb = string_tables(self.norb, self.k)
        self.tb = tb
        self.dim = dim = tb["dim"]
        self.K = tb["K"]
        dev = self.device

        # ---- backend
        self.backend = None
        if backend in ("auto", "numba"):
            try:
                self.backend = _NumbaBackend(dev, dtype)
            except Exception as e:  # noqa: BLE001
                if backend == "numba":
                    raise
                print(f"[gpu_energy] numba backend unavailable ({type(e).__name__}: {e}); using torch")
        if self.backend is None:
            self.backend = _TorchBackend(dev, dtype)

        with torch.cuda.device(dev):
            # ---- (norb, k) tables on the GPU (cached per device)
            g = _gpu_tables(self.norb, self.k, str(dev))
            self.g = g
            self.chunks = None
            if self.backend.name == "numba":
                self.pairs_b = g["pairs_nb"]
                if fuse and self.norb >= 4:                    # fused minor-index Givens passes
                    itc_ = 8 if dtype == torch.complex64 else 16
                    smem = int(smem_kb * 1024)
                    plan = plan_fusion(self.norb, self.k, itc_, smem)
                    t = fuse_t or (plan[0] if plan else None)
                    if t is not None and t < self.norb - 1:
                        self.chunks = _gpu_chunk_tables(self.norb, self.k, t, itc_, str(dev), smem, fuse_rows)
            else:
                self.pairs_b = g["pairs_t"]
            self.t_top = self.chunks["t"] if self.chunks is not None else 0

            self.eig_tol, self._max_mem_gb, self._tile_req, self._max_tile = eig_tol, max_mem_gb, tile, max_tile
            self._alloc_state = alloc_state
            self._load_hamiltonian(one_body, two_body, constant)

    def _load_hamiltonian(self, one_body, two_body, constant):
        """Factorize H, size the workspace, build the per-row coefficient table C and the tile buffers
        (and the state buffer on first use).  Depends on the molecule; everything else is reused."""
        dev, dim, dtype = self.device, self.dim, self.cdtype
        g = self.g
        with torch.cuda.device(dev):
            for a in ("C", "_Gbuf", "_Xbuf", "_out"):               # drop the previous molecule's buffers
                if hasattr(self, a):
                    delattr(self, a)
            torch.cuda.empty_cache()
            self.constant = float(constant)
            # ---- Hamiltonian factorization -> per-row coefficient tables C (dim, NL+1, K)
            gP, V, npos, E0, finfo = factorize_hamiltonian(one_body, two_body, self.k, self.eig_tol)
            self.fact_info = finfo
            self.E0 = E0
            self.nq, self.npos = V.shape[0], npos
            Vext = torch.as_tensor(np.concatenate([V, gP[None, :]], 0), device=dev, dtype=torch.float64)
            self.NLp = Vext.shape[0]

            # ---- memory plan
            itc = 8 if dtype == torch.complex64 else 16
            itr = itc // 2
            free, total = torch.cuda.mem_get_info(dev)
            psi_b = itc * dim * dim
            c_b = itr * dim * self.NLp * self.K
            have_psi = hasattr(self, "psi")
            if self._max_mem_gb is None:
                budget = free - 0.6e9 + (psi_b if have_psi else 0)
            else:
                budget = self._max_mem_gb * 1e9
            rest = budget - (psi_b if (self._alloc_state or have_psi) else 0) - c_b - 64e6
            self.mem_plan = dict(psi_gb=psi_b / 1e9, ctab_gb=c_b / 1e9, budget_gb=budget / 1e9,
                                 free_gb=free / 1e9)
            if rest < 0.15e9:
                raise MemoryError(f"norb={self.norb} k={self.k}: state {psi_b/1e9:.2f} GB + coefficient table "
                                  f"{c_b/1e9:.2f} GB exceed the budget {budget/1e9:.2f} GB")
            # energy tiles (preallocated): G (T^2 K) + X1, X2 (T^2 NL') complex -> 85% of the workspace.
            # (Diagonal-Coulomb phases need no workspace with the numba kernel; the GEMM fallback for a
            # non-diagonal J_ab sizes its row blocks from the free memory at call time.)
            per = itc * (self.K + 2 * self.NLp)
            frac = 0.85 if self.backend.name == "numba" else 0.55
            # cap: on the TITAN Xp the batched GEMM is 15-20% slower for T >= 1344 than for 640 <= T <= 1280
            T = self._tile_req or min(int(math.sqrt(frac * rest / per)), self._max_tile)
            T = max(32, min(dim, (T // 32) * 32 if T >= 64 else T))
            self.T = T
            self.mem_plan.update(tile=T, workspace_gb=rest / 1e9)

            # coefficient table, built in row chunks.  Diagonal entry of row a for operator kappa:
            # sum_p kappa_pp (n_p(a) - n_p(HF))  (HF shift, one spin), off-diagonal: kappa_P(a,e) sign(a,e)
            C = torch.empty((dim, self.NLp, self.K), device=dev, dtype=self.rdtype)
            exc_p, exc_s, occf = g["exc_p"], g["exc_s"], g["occf"]
            Vd = Vext[:, g["diagP"]]                                  # (NL', n)
            nhf = torch.zeros(self.norb, device=dev, dtype=torch.float64)
            nhf[: self.k] = 1.0
            step = max(1, int(2e8 // (self.NLp * self.K * 8)))
            for r0 in range(0, dim, step):
                r1 = min(dim, r0 + step)
                blk = Vext[:, exc_p[r0:r1]]                           # (NL', R, K-1)
                blk = blk.permute(1, 0, 2) * exc_s[r0:r1, None, :]
                C[r0:r1, :, :-1] = blk.to(self.rdtype)
                C[r0:r1, :, -1] = ((occf[r0:r1] - nhf) @ Vd.T).to(self.rdtype)
                del blk
            self.C = C
            del Vext
            if self._alloc_state and not have_psi:   # alloc_state=False: tables only (partial benchmarks)
                self.psi = torch.empty((dim, dim), device=dev, dtype=dtype)
                self.psi_b = self.backend.arr(self.psi)
            self._Gbuf = torch.empty(T * self.K * T, device=dev, dtype=dtype)
            self._Xbuf = torch.empty(2 * T * self.NLp * T, device=dev, dtype=dtype)
            self._out = torch.empty(((T + 31) // 32) ** 2 + 1, 4, device=dev, dtype=torch.float64)
            torch.cuda.synchronize(dev)

    def set_hamiltonian(self, one_body, two_body, constant, e_hf=None, e_ccsd=None, name=None):
        """Switch to another molecule with the same (norb, nelec), keeping the state buffer and every
        (norb, nelec) table (e.g. an RL loop over many norb-18 molecules on one 40 GB card)."""
        self.e_hf, self.e_ccsd, self.name = e_hf, e_ccsd, name
        self._load_hamiltonian(one_body, two_body, constant)
        return self

    def set_npz(self, path_or_name, ham_dir=None):
        """set_hamiltonian from rhf_hamiltonians/<name>.npz (norb and nelec must match)."""
        p = Path(path_or_name)
        if p.suffix != ".npz" or not p.exists():
            p = Path(ham_dir or _default_ham_dir()) / f"{path_or_name}.npz"
        d = np.load(p)
        if (int(d["norb"]), int(d["nelec_a"]), int(d["nelec_b"])) != (self.norb, self.k, self.k):
            raise ValueError(f"{p.stem}: (norb, nelec) differs from the engine's ({self.norb}, {self.nelec})")
        return self.set_hamiltonian(d["one_body"], d["two_body"], float(d["constant"]), e_hf=float(d["e_hf"]),
                                    e_ccsd=float(d["e_ccsd"]), name=p.stem)

    # ------------------------------------------------------------------ constructors
    @classmethod
    def from_npz(cls, path_or_name, ham_dir=None, **kw):
        p = Path(path_or_name)
        if p.suffix != ".npz" or not p.exists():
            p = Path(ham_dir or _default_ham_dir()) / f"{path_or_name}.npz"
        d = np.load(p)
        return cls(d["one_body"], d["two_body"], float(d["constant"]), int(d["norb"]),
                   (int(d["nelec_a"]), int(d["nelec_b"])), e_hf=float(d["e_hf"]),
                   e_ccsd=float(d["e_ccsd"]), name=p.stem, **kw)

    @classmethod
    def from_ffsim(cls, ham, norb, nelec, **kw):
        return cls(ham.one_body_tensor.real, ham.two_body_tensor.real, ham.constant, norb, nelec, **kw)

    def release(self):
        for a in ("psi", "psi_b", "C", "_Gbuf", "_Xbuf", "_out"):
            if hasattr(self, a):
                delattr(self, a)
        torch.cuda.empty_cache()

    # ------------------------------------------------------------------ parameters -> rotations
    def _prepare(self, U, Z, t1, connectivity):
        n = self.norb
        tonp = (lambda x: x.detach().cpu().numpy() if isinstance(x, torch.Tensor) else x)
        U = np.asarray(tonp(U), dtype=np.complex128)
        Z = np.asarray(tonp(Z), dtype=np.float64)
        t1 = None if t1 is None else np.asarray(tonp(t1))
        if U.ndim == 2:
            U = U[None]
        if Z.ndim == 2:
            Z = Z[None]
        n_reps = U.shape[0]
        assert U.shape == (n_reps, n, n) and Z.shape == (n_reps, n, n), (U.shape, Z.shape)
        if self.polar:
            W_, _, Vh = np.linalg.svd(U)
            U = W_ @ Vh
        maa, mab = interaction_masks(connectivity, n)
        if t1 is not None:
            from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
            F = orbital_rotation_from_t1_amplitudes(np.asarray(t1))
        else:
            F = np.eye(n)
        Ws = [U[0].conj().T] + [U[r].conj().T @ U[r - 1] for r in range(1, n_reps)] + [F @ U[-1]]
        dcs = []
        for r in range(n_reps):
            Jaa = Z[r] * maa
            Jab = Z[r] * mab
            # the spin-balanced UCJ needs a symmetric J_ab (ffsim validates the same with these tolerances);
            # psi = psi^T, which every step relies on, holds only then
            if not np.allclose(Jab, Jab.T, rtol=1e-5, atol=1e-8):
                raise ValueError("Z (J_ab part) must be symmetric for the spin-balanced UCJ operator")
            Jab = 0.5 * (Jab + Jab.T)
            dcs.append((np.triu(Jaa, 1), np.diag(Jaa).copy(), Jab))
        decs = []
        for W in Ws[1:]:
            app, D = givens_clustered(W, self.t_top)                 # application order, D applied first
            decs.append((split_segments(app, n, self.t_top), D))
        return Ws[0], dcs, decs

    # ------------------------------------------------------------------ state
    def _phase_factors(self, dc, D_next):
        """Per-string factor exp(i phi(s)) * prod_{p in s} D_next[p]  (complex128, GPU) and J_ab."""
        Jup, Jdiag, Jab = dc
        occf = self.g["occf"]
        Jup_t = torch.as_tensor(Jup, device=self.device)
        phi = ((occf @ Jup_t) * occf).sum(1) + 0.5 * occf @ torch.as_tensor(Jdiag, device=self.device)
        logd = torch.as_tensor(np.angle(D_next), device=self.device)
        ang = phi + occf @ logd
        mag = torch.exp(occf @ torch.as_tensor(np.log(np.abs(D_next)), device=self.device))
        return torch.polar(mag, ang), torch.as_tensor(Jab, device=self.device)

    def _phase_pass(self, rowfac, Jab, init: bool):
        """psi[a,b] = (init ? 1 : psi[a,b]) * rowfac[a] rowfac[b] exp(i occ_a^T Jab occ_b)."""
        occf = self.g["occf"]
        diagonal_only = not bool(torch.count_nonzero(Jab - torch.diag(torch.diagonal(Jab))))
        if diagonal_only and self.backend.name == "numba":
            ez = torch.polar(torch.ones(self.norb, device=self.device, dtype=torch.float64), torch.diagonal(Jab))
            self.backend.phase_diag(self.psi_b, self.g["strs_nb"], rowfac.contiguous(), ez, self.norb, init)
            return
        JO = Jab @ occf.T                                          # (n, dim) fp64
        dim = self.dim
        itc = 8 if self.cdtype == torch.complex64 else 16
        free = torch.cuda.mem_get_info(self.device)[0] - 0.3e9     # theta+ones+polar+cast per element
        R = max(1, min(dim, int(0.8 * free / (dim * (32 + itc)))))
        for r0 in range(0, dim, R):
            r1 = min(dim, r0 + R)
            th = occf[r0:r1] @ JO if not diagonal_only else (occf[r0:r1] * torch.diagonal(Jab)) @ occf.T
            f = torch.polar(torch.ones_like(th), th)
            f.mul_(rowfac[r0:r1, None])
            f.mul_(rowfac[None, :])
            if init:
                self.psi[r0:r1].copy_(f)
            else:
                self.psi[r0:r1].mul_(f.to(self.cdtype))
            del th, f

    def _orbital_rotation(self, segs):
        """psi -> A psi A^T with A = Lambda^k(W) (segments in application order; D of W already applied)."""
        be = self.backend
        if self.chunks is not None:
            # X <- X A^T (columns), X <- X^T (= A X since X was symmetric), X <- X A^T  =>  A X A^T
            useg = be.upload_segments(segs)
            be.givens_minor(self.psi_b, self.dim, useg, self.pairs_b, self.chunks)
            be.transpose(self.psi_b, self.dim)
            be.givens_minor(self.psi_b, self.dim, useg, self.pairs_b, self.chunks)
            return
        # unfused: one pass per rotation on the rows: A X, transpose (= X^T A^T), A X^T A^T = A X A^T
        rots = [r for _, seg in segs for r in seg]
        be.givens_sequence(self.psi_b, self.dim, rots, self.pairs_b)
        be.transpose(self.psi_b, self.dim)
        be.givens_sequence(self.psi_b, self.dim, rots, self.pairs_b)

    def _build_state(self, U, Z, t1, connectivity):
        tm = {}
        t0 = time.perf_counter()
        W1, dcs, decs = self._prepare(U, Z, t1, connectivity)
        tm["prepare"] = time.perf_counter() - t0
        sync = torch.cuda.synchronize
        # Slater determinant from HF
        t0 = time.perf_counter()
        occl = self.g["occ_list"]
        W1t = torch.as_tensor(W1, device=self.device)
        v = torch.linalg.det(W1t[occl][:, :, : self.k])           # (dim,) complex128
        rowfac, Jab = self._phase_factors(dcs[0], decs[0][1])
        self._phase_pass(v * rowfac, Jab, init=True)
        sync(self.device)
        tm["init+phase1"] = time.perf_counter() - t0
        tm["rotations"] = 0.0
        tm["phases"] = 0.0
        n_reps = len(dcs)
        for r in range(n_reps):
            t0 = time.perf_counter()
            self._orbital_rotation(decs[r][0])
            sync(self.device)
            tm["rotations"] += time.perf_counter() - t0
            if r + 1 < n_reps:
                t0 = time.perf_counter()
                rowfac, Jab = self._phase_factors(dcs[r + 1], decs[r + 1][1])
                self._phase_pass(rowfac, Jab, init=False)
                sync(self.device)
                tm["phases"] += time.perf_counter() - t0
        return tm

    # ------------------------------------------------------------------ energy
    def _x_tile(self, a0, a1, b0, b1, slot):
        """X[A, :, B] (nA, NL', nB) complex view into the workspace."""
        nA, nB = a1 - a0, b1 - b0
        K, NLp = self.K, self.NLp
        idx = self.g["src"][a0:a1].reshape(-1)
        G = self._Gbuf[: nA * K * nB].view(nA * K, nB)
        torch.index_select(self.psi[:, b0:b1], 0, idx, out=G)
        Gr = torch.view_as_real(G).view(nA, K, 2 * nB)
        off = slot * self.T * NLp * self.T
        Xc = self._Xbuf[off: off + nA * NLp * nB]
        Xr = torch.view_as_real(Xc).view(nA, NLp, 2 * nB)
        torch.bmm(self.C[a0:a1], Gr, out=Xr)
        return Xc.view(nA, NLp, nB)

    def _energy_of_psi(self):
        dim, T = self.dim, self.T
        # the exact state is symmetric; remove the antisymmetric rounding part, which would otherwise
        # enter the symmetric-psi energy formula at first order
        self.backend.symmetrize(self.psi_b, dim)
        tiles = [(i, min(dim, i + T)) for i in range(0, dim, T)]
        acc = torch.zeros(4, device=self.device, dtype=torch.float64)    # e2 (x0.5 later), e2neg, e1, norm
        prev_tf32 = torch.backends.cuda.matmul.allow_tf32
        torch.backends.cuda.matmul.allow_tf32 = False
        try:
            for i, (a0, a1) in enumerate(tiles):
                for j in range(i, len(tiles)):
                    b0, b1 = tiles[j]
                    X1 = self._x_tile(a0, a1, b0, b1, 0)
                    diag = i == j
                    X2 = X1 if diag else self._x_tile(b0, b1, a0, a1, 1)
                    P1 = self.psi[a0:a1, b0:b1]
                    P2 = P1 if diag else self.psi[b0:b1, a0:a1]
                    nb = self.backend.epilogue(X1, X2, P1, P2, self.nq, self.npos, diag, self._out)
                    s = self._out[:nb].sum(0)
                    f = 1.0 if diag else 2.0
                    acc[0] += f * (s[0] - s[1])
                    acc[2] += s[2]
                    acc[3] += s[3]
        finally:
            torch.backends.cuda.matmul.allow_tf32 = prev_tf32
        e2, _, e1, nrm = acc.tolist()
        E1 = 2.0 * e1
        E2 = 0.5 * e2
        self.last = dict(E1=E1, E2=E2, norm=nrm, E0=self.E0)
        return self.constant + self.E0 + (E1 + E2) / nrm

    # ------------------------------------------------------------------ public API
    @torch.no_grad()
    def energy(self, U, Z, t1=None, connectivity: str = "square") -> float:
        """Exact LUCJ energy (Hartree) of parameters U (n_reps,n,n) complex, Z (n_reps,n,n) real."""
        with torch.cuda.device(self.device):
            t0 = time.perf_counter()
            tm = self._build_state(U, Z, t1, connectivity)
            t1_ = time.perf_counter()
            E = self._energy_of_psi()
            tm["energy"] = time.perf_counter() - t1_
            tm["total"] = time.perf_counter() - t0
            self.timing = tm
            return float(E)

    def energies(self, params, t1=None, connectivity: str = "square") -> list[float]:
        """params: iterable of (U, Z) or (U, Z, t1); t1 given here is used when an item has none."""
        out = []
        for item in params:
            if len(item) == 3:
                U, Z, tt = item
            else:
                (U, Z), tt = item, None
            out.append(self.energy(U, Z, t1=tt if tt is not None else t1, connectivity=connectivity))
        return out

    def corr_frac(self, E):
        return corr_frac(E, self.e_hf, self.e_ccsd)

    @torch.no_grad()
    def state(self, U, Z, t1=None, connectivity: str = "square") -> torch.Tensor:
        """The (dim, dim) CI matrix (internal buffer; clone it before the next call)."""
        with torch.cuda.device(self.device):
            self._build_state(U, Z, t1, connectivity)
            return self.psi

    @torch.no_grad()
    def energy_of_state(self, psi) -> float:
        """Energy of a given symmetric CI matrix / ffsim vector (copied into the internal buffer)."""
        with torch.cuda.device(self.device):
            psi = torch.as_tensor(psi).reshape(self.dim, self.dim)
            self.psi.copy_(psi)
            return float(self._energy_of_psi())


@lru_cache(maxsize=None)
def _gpu_chunk_tables(norb: int, k: int, t: int, itemsize: int, device: str,
                      smem_bytes: int = 48 * 1024, rows: int | None = None) -> dict:
    """Device tables of the fused minor-index kernel; R rows of the largest chunk per block."""
    from numba import cuda
    ct = chunk_tables(norb, k, t)
    dev = torch.device(device)
    keep = {}
    stride = int(ct["max_len"])
    pair_b = 4 * ct["max_kl_len"] + 16
    R = rows or max(1, min(16, (smem_bytes - pair_b) // (stride * itemsize)))
    u_off = R * stride * (itemsize // 2)                 # offset of the uint16 region, in uint16 units
    out = dict(t=t, nch=len(ct["starts"]), max_len=ct["max_len"], R=int(R), stride=stride, u_off=int(u_off),
               smem_minor=int(R * stride * itemsize + pair_b))
    for key in ("starts", "lens", "kls", "cnt", "lo16", "hi16", "kl_base", "kl_len", "off_rel"):
        keep[key] = torch.as_tensor(np.ascontiguousarray(ct[key]).view(np.int16) if ct[key].dtype == np.uint16
                                    else np.ascontiguousarray(ct[key]), device=dev)
        arr = cuda.as_cuda_array(keep[key])
        out[key] = arr.view(np.uint16) if ct[key].dtype == np.uint16 else arr
    out["_keep"] = keep                                  # keep the torch storage alive
    return out


@lru_cache(maxsize=None)
def _gpu_tables(norb: int, k: int, device: str) -> dict:
    tb = string_tables(norb, k)
    dev = torch.device(device)
    occf = torch.as_tensor(tb["occ"], device=dev, dtype=torch.float64)
    pairs_t = [(torch.as_tensor(lo, device=dev, dtype=torch.int64),
                torch.as_tensor(hi, device=dev, dtype=torch.int64)) for lo, hi in tb["pairs"]]
    pairs_i32 = [(torch.as_tensor(lo, device=dev), torch.as_tensor(hi, device=dev)) for lo, hi in tb["pairs"]]
    g = dict(occf=occf, occ_list=torch.as_tensor(tb["occ_list"], device=dev),
             strs=torch.as_tensor(tb["strs"], device=dev),
             src=torch.as_tensor(tb["src"], device=dev),
             exc_p=torch.as_tensor(tb["exc_p"], device=dev),
             exc_s=torch.as_tensor(tb["exc_s"], device=dev, dtype=torch.float64),
             diagP=torch.as_tensor(tb["diagP"], device=dev), pairs_t=pairs_t, pairs_i32=pairs_i32)
    try:
        from numba import cuda
        g["pairs_nb"] = [(cuda.as_cuda_array(lo), cuda.as_cuda_array(hi)) for lo, hi in pairs_i32]
        g["strs_nb"] = cuda.as_cuda_array(g["strs"])
    except Exception:  # noqa: BLE001
        g["pairs_nb"] = None
    return g


# ======================================================================= CLI: GPU drop-in for the CPU energy job

def main():
    """Exact energies of every task in an energy-task pickle ({"tasks": [((name, cand), name, U, Z, t1), ...]}),
    one engine per molecule; writes per_molecule[name][cand] = {E, corr_frac, t} like the CPU energy job.

      CUDA_VISIBLE_DEVICES=0 python -m pretrain.rl.gpu_energy --tasks runs_ot/energy_tasks/n17_rl1.pkl \\
          --out pretrain/opt_true/results/energy_n17_rl1_gpu.json [--dtype c8|c16] [--connectivity square]
    """
    import argparse
    import json
    import pickle

    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ham-dir", default=str(_default_ham_dir()))
    ap.add_argument("--dtype", choices=["c8", "c16"], default="c8")
    ap.add_argument("--connectivity", default="square")
    ap.add_argument("--max-mem-gb", type=float, default=None)
    args = ap.parse_args()
    dt = torch.complex64 if args.dtype == "c8" else torch.complex128
    tasks = pickle.load(open(args.tasks, "rb"))["tasks"]
    by_mol: dict = {}
    for key, name, U, Z, t1 in tasks:
        by_mol.setdefault(name, []).append((key, U, Z, t1))
    res: dict = {}
    t_all = time.time()
    for name, items in by_mol.items():
        eng = LUCJEnergyGPU.from_npz(name, ham_dir=args.ham_dir, dtype=dt, max_mem_gb=args.max_mem_gb)
        for key, U, Z, t1 in items:
            E = eng.energy(U, Z, t1=t1, connectivity=args.connectivity)
            cf = eng.corr_frac(E)
            cand = key[1] if isinstance(key, tuple) else str(key)
            res.setdefault(name, {})[cand] = {"E": E, "corr_frac": cf, "t": eng.timing["total"]}
            print(f"  {name:22s} {cand:28s} corr% {100 * cf:8.3f}  E {E:.10f}  ({eng.timing['total']:.1f}s)",
                  flush=True)
        eng.release()
        del eng
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump({"per_molecule": res, "args": vars(args), "engine": "pretrain.rl.gpu_energy",
               "gpu": torch.cuda.get_device_name()}, open(args.out, "w"), indent=1)
    print(f"-> {args.out}  ({time.time() - t_all:.0f}s)")


if __name__ == "__main__":
    main()
