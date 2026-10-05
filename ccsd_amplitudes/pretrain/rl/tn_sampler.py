"""Perfect sampling of configurations from the dm-engine MPS (pretrain/rl/tn_energy.py) -- TN-sampled SQD / QSCI
beyond exact state vectors (norb 29 = 58 qubits and up).

Public API
    ev = LUCJEnergyTN(one_body, two_body, constant, norb, nelec, max_bond=256, device="cuda", name=name)
    mps, Fr, info = build_mps(ev, U, Z, t1)                 # the reward engine's MPS, built once
    E, einfo = mps_energy(ev, mps, Fr)                      # optional: block2 <H> of that MPS (== ev.energy())
    smp = MPSSampler(mps)                                   # right-canonicalizes (block QR), GPU tensors
    occ, logp = smp.sample(100_000, seed=0)                 # occ (N, norb) uint8, site state s = n_alpha + 2 n_beta
    ints = occ_to_ints(occ)                                 # ffsim / qiskit integer bitstrings (alpha = low bits)
    h, g = sampled_basis_integrals(one_body, two_body, Fr)  # the Hamiltonian in the basis of the bitstrings
    batches = sqd_batches(ints, norb, nelec, seed=0)        # Lin et al.'s SQD subsampling (qiskit-addon-sqd 0.13.1)
    E_b = [solve_batch(h, g, constant, b, norb, nelec) for b in batches]   # QSCI energy per 4,000-sample batch

Which orbital basis the bitstrings refer to (precisely)
  * MO basis: the active-space RHF MOs phi_p of rhf_hamiltonians/<name>.npz (one_body h, two_body (pq|rs), chemists'
    notation; AO coefficients = rhf_dataset/<name>.npz "mo_coeff", the active MOs, signs aligned with the dataset).
  * The LUCJ state of the energy jobs is  |psi> = O(F U_L) e^{iJ_L} ... O(U_2^dag U_1) e^{iJ_1} O(U_1^dag) |HF>,
    F = orbital_rotation_from_t1_amplitudes(t1) (ffsim's make_ucj_op(Z, U, "square", t1)), with O(W) the orbital
    rotation a^dag_i -> sum_j W_ji a^dag_j (ffsim convention).
  * The engine works in the split-localized basis S (n x n, real orthogonal, block-diagonal w.r.t. occupied/virtual
    MOs; column j = site j, Fiedler order; tn_energy.split_localized_basis).  Its MPS holds the coefficients
        c(x) = <x; chi | psi>,   chi_j = sum_p (Fr)_pj phi_p,   Fr = F @ S   (returned by ev.state / build_mps),
    i.e. c = O(Fr)^dag |psi> = ffsim.apply_orbital_rotation(psi, Fr^T, norb, nelec).  The final t1 rotation F and
    the basis change S are NOT applied to the MPS; they are folded into the Hamiltonian, which the engine evaluates
    as O(Fr)^dag H O(Fr).  A sampled bitstring therefore gives the occupations n_{j,sigma} of the orbitals chi_j
    (AO coefficients C_mo_active @ F @ S), NOT of the MOs.  Fr is real (t1 and S real), so the matching integrals
        h^chi = Fr^T h Fr,   (jk|lm)^chi = sum Fr_pj Fr_qk Fr_rl Fr_sm (pq|rs),   constant unchanged
    are real (sampled_basis_integrals == tn_energy.rotate_hamiltonian) and SQD must be run with them.  Bit j of the
    alpha (beta) integer = n_{j,alpha} (n_{j,beta}) = site j of the MPS.
  * SQD/QSCI energies depend on the basis (the subspace is a product of sampled alpha x beta strings).  Lin et al.
    sample the circuit's qubits, i.e. the MO basis (their circuit includes the final orbital rotation).  chi-basis
    QSCI is an equally variational estimator but a different one; MO-basis sampling from this MPS would require
    applying O(Fr) to it, a volume-law rotation (see the tn_energy docstring).  docs/followups/sampler.md quantifies
    the MO-vs-chi difference at norb 15-16 with exact samples.
  * Amplitude convention.  The MPS is in the interleaved Jordan-Wigner order (0a, 0b, 1a, 1b, ...; creation operators
    ascending left to right), ffsim/pyscf vectors in the alpha-block-then-beta order.  For the same occupations
        c_ffsim(x) = phase * (-1)^{sum_p n_beta(p) sum_{q>p} n_alpha(q)} * c_mps(x)
    with one global phase (verified on random systems, pretrain/rl/tests/test_tn_sampler.py).  Probabilities are equal.

Sampling (perfect / direct sampling, Ferris & Vidal PRB 85, 165146 (2012)): with the MPS right-canonical (centre at
site 0, block QR of the engine), p(s_i | s_<i) = |L_i A_i[s_i]|^2 / |L_i|^2 exactly; the sampler draws site by site
from left to right for a batch of samples at once (dense GPU matmuls (B, chi) @ (chi, 4 chi); the U(1)xU(1) labels
make the rows of L_i zero outside one sector, which costs flops but needs no bookkeeping).  Each sample is an exact
draw from |c|^2 / <c|c> of the truncated MPS (up to the complex64 isometry error of the canonical form, ~1e-6
relative); `logp` is its log-probability.  The configurations conserve N_alpha and N_beta exactly.

Exact references (validation at norb <= 18): LUCJStateGPU is pretrain.rl.gpu_energy.LUCJEnergyGPU with the final
orbital rotation followed by a basis change B (columns in MO coordinates): its CI matrix is O(B)^dag |psi>, built by
folding B^dag into the last Givens network (no extra pass).  With B = Fr it is exactly the state the MPS
approximates; energy() with the chi-basis integrals reproduces the MO-basis energy (checked).
sample_ci_matrix() draws exact samples from a CI matrix on the GPU.
"""
from __future__ import annotations

import math
import time

import numpy as np
import torch

from pretrain.rl.gpu_energy import LUCJEnergyGPU, givens_clustered, split_segments

SQD_VERSION_TESTED = "0.13.1"

# SQD settings of Lin et al. (wan-hsuan-lucj scripts/quimb/*/lucj_compressed_t2.py -> quimb_task/
# lucj_sqd_quimb_task_sci.py): 10^5 samples, max_dim = samples_per_batch = 4000, 10 batches, max_iterations 1
# (no configuration recovery: the noiseless samples are postselected only), symmetrize_spin, solve_sci_batch.
LIN_SQD = dict(samples_per_batch=4000, num_batches=10, energy_tol=1e-5, occupancies_tol=1e-3,
               carryover_threshold=1e-3, max_iterations=1, symmetrize_spin=True, max_dim=4000)
LIN_SHOTS = 100_000


# ------------------------------------------------------------------------------------------- bitstring helpers

def occ_to_strings(occ: np.ndarray):
    """Site states (N, n) in {0,1,2,3} (s = n_alpha + 2 n_beta) -> (alpha int64 (N,), beta int64 (N,)),
    bit j = occupation of orbital/site j."""
    occ = np.asarray(occ)
    n = occ.shape[1]
    w = np.left_shift(np.int64(1), np.arange(n, dtype=np.int64))
    a = ((occ & 1).astype(np.int64) * w).sum(1)
    b = ((occ >> 1).astype(np.int64) * w).sum(1)
    return a, b


def strings_to_occ(a: np.ndarray, b: np.ndarray, n: int) -> np.ndarray:
    a = np.asarray(a, dtype=np.int64)[:, None]
    b = np.asarray(b, dtype=np.int64)[:, None]
    j = np.arange(n, dtype=np.int64)[None, :]
    return (((a >> j) & 1) + 2 * ((b >> j) & 1)).astype(np.uint8)


def occ_to_ints(occ: np.ndarray) -> np.ndarray:
    """Site states -> one integer per sample in the ffsim BitstringType.INT / qiskit BitArray convention:
    bits 0..n-1 = alpha orbitals, bits n..2n-1 = beta orbitals (the alpha part is the right half of the qiskit
    bitstring [b_{n-1} .. b_0 a_{n-1} .. a_0], as qiskit_addon_sqd expects).  n <= 31."""
    occ = np.asarray(occ)
    n = occ.shape[1]
    if 2 * n > 63:
        raise ValueError("2 * norb > 63 bits: use occ_to_strings")
    a, b = occ_to_strings(occ)
    return a | (b << np.int64(n))


def ints_to_strings(ints: np.ndarray, n: int):
    ints = np.asarray(ints, dtype=np.int64)
    mask = np.int64((1 << n) - 1)
    return ints & mask, (ints >> np.int64(n)) & mask


def ffsim_sign(occ: np.ndarray) -> np.ndarray:
    """(-1)^{sum_p n_beta(p) sum_{q>p} n_alpha(q)}: interleaved (MPS) -> alpha-block (ffsim) amplitude sign."""
    occ = np.asarray(occ).astype(np.int64)
    na, nb = occ & 1, occ >> 1
    suf = np.cumsum(na[:, ::-1], axis=1)[:, ::-1] - na
    return 1 - 2 * ((nb * suf).sum(1) & 1)


# --------------------------------------------------------------------------------------------------- basis

def sampled_basis_integrals(one_body, two_body, Fr):
    """Integrals of O(Fr)^dag H O(Fr): the Hamiltonian in the orbital basis chi_j = sum_p Fr_pj phi_p of the
    sampled bitstrings.  Real for real Fr (always the case here: Fr = final(t1) @ S).  The constant is unchanged."""
    from pretrain.rl.tn_energy import rotate_hamiltonian
    h, g = rotate_hamiltonian(one_body, two_body, Fr)
    if np.isrealobj(Fr) or (np.abs(np.imag(Fr)).max() < 1e-14):
        assert np.abs(h.imag).max() < 1e-10 and np.abs(g.imag).max() < 1e-10
        return np.ascontiguousarray(h.real), np.ascontiguousarray(g.real)
    return h, g


def sampled_basis(ev, t1):
    """Fr = final(t1) @ S of a tn_energy evaluator (the columns are the sampled orbitals in MO coordinates)."""
    if t1 is None:
        return np.asarray(ev.S)
    from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
    return orbital_rotation_from_t1_amplitudes(np.asarray(t1, dtype=np.float64)) @ np.asarray(ev.S)


# ----------------------------------------------------------------------------------------- MPS: build / energy

def build_mps(ev, U, Z, t1=None, max_bond: int | None = None):
    """The engine's MPS (ev.state) plus its truncation record.  Returns (mps, Fr, info)."""
    t0 = time.time()
    chi = ev.max_bond if max_bond is None else int(max_bond)
    mps, Fr = ev.state(U, Z, t1, chi)
    if mps.device.type == "cuda":
        torch.cuda.synchronize(mps.device)
    info = {"chi": chi, "method": ev.method, "discarded_sum": float(mps.discarded), "discarded_max": float(mps.max_disc),
            "n_trunc": int(mps.n_trunc), "discarded_zip_sum": float(getattr(mps, "discarded_zip", 0.0)),
            "max_bond": int(max(mps.bond_dims())), "bond_dims": [int(x) for x in mps.bond_dims()],
            "t_state": time.time() - t0, "settings": ev.settings()}
    return mps, Fr, info


_MPO_CACHE: dict = {}


def mps_energy(ev, mps, Fr):
    """<mps|H|mps>/<mps|mps> with the engine's block2 path (the same numbers as ev.energy, without rebuilding the
    state).  The MPO is cached per (evaluator, Fr)."""
    from pretrain.rl.tn_energy import rotate_hamiltonian
    t0 = time.time()
    mps.move_center(0)
    nrm2 = float((abs(mps.A[0]) ** 2).sum())
    key = (id(ev), hash(np.asarray(Fr).tobytes()))
    if key not in _MPO_CACHE:
        _MPO_CACHE.clear()
        h, g = rotate_hamiltonian(ev.h, ev.g, Fr)
        _MPO_CACHE[key] = ev.b2.mpo(h, g, ev.const)
    t1 = time.time()
    bm = ev.b2.to_block2(mps)
    e = ev.b2.expectation(bm, _MPO_CACHE[key]) / nrm2
    del bm
    return float(e.real), {"imag": float(e.imag), "t_mpo": t1 - t0, "t_expect": time.time() - t1}


# --------------------------------------------------------------------------------------------------- sampler

class MPSSampler:
    """Direct sampling and amplitude evaluation for a tn_energy MPS (GpuSymMPS, NpSymMPS or SymMPS).

    The constructor moves the orthogonality centre to site 0 with the engine's block QR (sites 1..n-1 become right
    isometries; the MPS object is modified in place, its state is unchanged) and copies the tensors to `device`
    (default: the MPS device, or CUDA if available for numpy engines) in `dtype` (default: the MPS precision).

        smp = MPSSampler(mps)
        occ, logp = smp.sample(n, seed=0)          # occ (n, norb) uint8, logp = log p_mps(x) (float64)
        amp = smp.amplitudes(occ)                  # complex128 amplitudes of the normalized MPS (interleaved order)
        p = smp.probabilities(occ)
    """

    def __init__(self, mps, device=None, dtype=None, canonicalize: bool = True, check: bool = True):
        if canonicalize:
            # full block-QR sweep right and back: every site 1..n-1 becomes an exact right isometry even if the
            # engine left zero-weight (non-orthonormal) bond directions behind (zip-up engine with cutoff 0)
            mps.move_center(mps.n - 1)
            mps.move_center(0)
        self.n, self.nelec = int(mps.n), tuple(int(x) for x in mps.nelec)
        A0 = mps.A[0]
        if device is None:
            if isinstance(A0, torch.Tensor):
                device = A0.device
            else:
                device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device)
        if dtype is None:
            dtype = A0.dtype if isinstance(A0, torch.Tensor) else \
                (torch.complex64 if A0.dtype == np.complex64 else torch.complex128)
        self.dtype = dtype
        self.A = [torch.as_tensor(a).to(self.device, dtype).contiguous() for a in mps.A]
        self.q = [np.asarray(q.cpu().numpy() if isinstance(q, torch.Tensor) else q) for q in mps.q]
        self.norm2 = float((self.A[0].abs() ** 2).sum())
        self.bond_dims = [int(a.shape[2]) for a in self.A[:-1]]
        self.isometry_err = self._isometry_error() if (check and canonicalize) else None

    @torch.no_grad()
    def _isometry_error(self) -> float:
        err = 0.0
        for A in self.A[1:]:
            Dl = A.shape[0]
            M = A.reshape(Dl, -1)
            G = M @ M.mH
            err = max(err, float((G - torch.eye(Dl, dtype=G.dtype, device=G.device)).abs().max()))
        return err

    # -------------------------------------------------------------- sampling
    @torch.no_grad()
    def sample(self, n_samples: int, seed: int = 0, batch: int | None = None, return_logp: bool = True):
        """n_samples exact draws from |c(x)|^2 / <c|c>.  Returns (occ (N, n) uint8, logp (N,) float64 or None).
        Uniforms come from numpy's PCG64 (np.random.default_rng(seed)) on the host, one (N, n) block, so the draws
        are identical on every device and for every batch size.  (torch's CUDA generator was NOT usable: with
        torch 2.4, consecutive float64 torch.rand calls on CUDA repeated parts of the stream -- 24,576 duplicates
        in 4 x (65,536 x 6) draws -- and biased a chi-square test at 5 sigma.)"""
        n = self.n
        if batch is None:
            chi = max(self.bond_dims) if self.bond_dims else 1
            batch = int(max(1024, min(1 << 16, (1 << 27) // (4 * chi))))     # V (B, 4 chi) <= 2^27 elements
        rng = np.random.default_rng(seed)
        occ_out = np.empty((n_samples, n), dtype=np.uint8)
        lp_out = np.empty(n_samples, dtype=np.float64) if return_logp else None
        for b0 in range(0, n_samples, batch):
            B = min(batch, n_samples - b0)
            u = torch.as_tensor(rng.random((B, n)), device=self.device)
            L = torch.ones((B, 1), dtype=self.dtype, device=self.device)
            lp = torch.zeros(B, dtype=torch.float64, device=self.device)
            occ = torch.empty((B, n), dtype=torch.uint8, device=self.device)
            ar = torch.arange(B, device=self.device)
            for i in range(n):
                A = self.A[i]
                Dl, _, Dr = A.shape
                V = (L @ A.reshape(Dl, 4 * Dr)).reshape(B, 4, Dr)
                w = (V.real.square() + V.imag.square()).sum(-1).to(torch.float64)      # (B, 4)
                c = torch.cumsum(w, dim=1)
                r = u[:, i] * c[:, 3]
                s = (r[:, None] >= c[:, :3]).sum(1)
                ws = w[ar, s]
                bad = ws <= 0                                                        # rounding at a boundary
                if bool(bad.any()):
                    s = torch.where(bad, w.argmax(1), s)
                    ws = w[ar, s]
                lp += torch.log(ws) - torch.log(c[:, 3])
                L = V[ar, s] / torch.sqrt(ws).to(V.real.dtype)[:, None]
                occ[:, i] = s.to(torch.uint8)
            occ_out[b0:b0 + B] = occ.cpu().numpy()
            if return_logp:
                lp_out[b0:b0 + B] = lp.cpu().numpy()
        return occ_out, lp_out

    # -------------------------------------------------------------- amplitudes
    @torch.no_grad()
    def amplitudes(self, occ: np.ndarray, batch: int = 1 << 15) -> np.ndarray:
        """Amplitudes c_mps(x) / sqrt(<c|c>) (complex128) of given site-state rows, in the MPS (interleaved JW)
        convention; multiply by ffsim_sign(occ) (and one global phase) for ffsim's convention.  Zero for
        configurations outside the (N_alpha, N_beta) sector."""
        occ = np.asarray(occ, dtype=np.uint8)
        N, n = occ.shape
        assert n == self.n
        out = np.empty(N, dtype=np.complex128)
        for b0 in range(0, N, batch):
            o = torch.as_tensor(occ[b0:b0 + batch].astype(np.int64), device=self.device)
            B = o.shape[0]
            L = torch.ones((B, 1), dtype=self.dtype, device=self.device)
            logn = torch.zeros(B, dtype=torch.float64, device=self.device)
            dead = torch.zeros(B, dtype=torch.bool, device=self.device)
            for i in range(n):
                A = self.A[i]
                Ln = torch.zeros((B, A.shape[2]), dtype=self.dtype, device=self.device)
                for sv in range(4):
                    m = (o[:, i] == sv).nonzero().squeeze(1)
                    if m.numel():
                        Ln[m] = L[m] @ A[:, sv, :]
                nr = torch.linalg.vector_norm(Ln, dim=1).to(torch.float64)
                z = nr <= 0
                dead |= z
                nr = torch.where(z, torch.ones_like(nr), nr)
                L = Ln / nr.to(Ln.real.dtype)[:, None]
                logn += torch.log(nr)
            amp = L[:, 0].to(torch.complex128) * torch.exp(logn - 0.5 * math.log(self.norm2)).to(torch.complex128)
            amp[dead] = 0.0
            out[b0:b0 + B] = amp.cpu().numpy()
        return out

    def probabilities(self, occ: np.ndarray, batch: int = 1 << 15) -> np.ndarray:
        return np.abs(self.amplitudes(occ, batch)) ** 2


# --------------------------------------------------------------------- exact references (GPU, norb <= 18-19)

class LUCJStateGPU(LUCJEnergyGPU):
    """pretrain.rl.gpu_energy.LUCJEnergyGPU whose final orbital rotation is followed by the basis change B (n x n
    unitary, columns = new orbitals in MO coordinates): the CI matrix is O(B)^dag |psi>, i.e. the state expressed in
    the orbitals chi_j = sum_p B_pj phi_p.  B is folded into the last Givens network: W_final = B^dag F U_L (no extra
    pass).  Pass the integrals in the same basis (sampled_basis_integrals) to get energies.  basis=None: the plain
    MO-basis engine."""

    def __init__(self, one_body, two_body, constant, norb, nelec, basis=None, **kw):
        self.basis = None if basis is None else np.asarray(basis, dtype=np.complex128)
        super().__init__(one_body, two_body, constant, norb, nelec, **kw)

    def _prepare(self, U, Z, t1, connectivity):
        W1, dcs, decs = super()._prepare(U, Z, t1, connectivity)
        if self.basis is None:
            return W1, dcs, decs
        tonp = (lambda x: x.detach().cpu().numpy() if isinstance(x, torch.Tensor) else x)
        Uk = np.asarray(tonp(U), dtype=np.complex128)
        Uk = Uk[None] if Uk.ndim == 2 else Uk
        UL = Uk[-1]
        if self.polar:
            W_, _, Vh = np.linalg.svd(UL)
            UL = W_ @ Vh
        if t1 is not None:
            from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
            F = orbital_rotation_from_t1_amplitudes(np.asarray(tonp(t1)))
        else:
            F = np.eye(self.norb)
        app, D = givens_clustered(self.basis.conj().T @ F @ UL, self.t_top)
        decs = list(decs)
        decs[-1] = (split_segments(app, self.norb, self.t_top), D)
        return W1, dcs, decs


@torch.no_grad()
def sample_ci_matrix(psi: torch.Tensor, n_samples: int, seed: int = 0):
    """Exact samples (ia, ib) (row = alpha string address, column = beta, ffsim/pyscf order) from |psi|^2 of a CI
    matrix on its device: inverse-CDF sampling on the float64 cumulative sum (torch.multinomial is limited to 2^24
    categories); uniforms from numpy PCG64 on the host (see MPSSampler.sample).  Returns numpy int64 arrays."""
    dim_a, dim_b = psi.shape
    p = psi.reshape(-1)
    c = p.real.to(torch.float64).square()
    c += p.imag.to(torch.float64).square()
    c = torch.cumsum(c, 0)
    u = torch.as_tensor(np.random.default_rng(seed).random(n_samples), device=psi.device) * c[-1]
    idx = torch.searchsorted(c, u, right=True).clamp_(max=c.numel() - 1)
    del c
    idx = idx.cpu().numpy()
    return idx // dim_b, idx % dim_b


def ci_strings(norb: int, k: int) -> np.ndarray:
    from pyscf.fci import cistring
    return np.asarray(cistring.make_strings(range(norb), k), dtype=np.int64)


def ci_addresses(norb: int, k: int, strs: np.ndarray) -> np.ndarray:
    from pyscf.fci import cistring
    return np.asarray(cistring.strs2addr(norb, k, np.asarray(strs, dtype=np.int64)), dtype=np.int64)


# ---------------------------------------------------------------------------------------- SQD (Lin et al.)

def bit_array_from_ints(ints: np.ndarray, norb: int):
    """qiskit BitArray of integer bitstrings (alpha = low bits), via counts."""
    from qiskit.primitives import BitArray
    vals, cnt = np.unique(np.asarray(ints, dtype=np.int64), return_counts=True)
    return BitArray.from_counts({int(v): int(c) for v, c in zip(vals, cnt)}, num_bits=2 * norb)


def sqd_batches(ints: np.ndarray, norb: int, nelec, seed=0, samples_per_batch: int = LIN_SQD["samples_per_batch"],
                num_batches: int = LIN_SQD["num_batches"], max_dim: int | None = LIN_SQD["max_dim"],
                symmetrize_spin: bool = LIN_SQD["symmetrize_spin"]):
    """The CI-string pairs (strs_a, strs_b) of the first SQD iteration exactly as
    qiskit_addon_sqd.fermion.diagonalize_fermionic_hamiltonian builds them (postselection on the Hamming weights,
    subsample(samples_per_batch, num_batches) without replacement per batch, spin symmetrization, max_dim
    truncation by marginal counts).  With max_iterations = 1 (Lin et al.) the SQD energy is
    min_b solve_batch(..., batches[b]).  Calls the library's own (private, 0.13.1) _prepare_ci_strings, so the
    subspaces are identical to the library's for the same seed; the subspace dimension of batch b is
    len(strs_a) * len(strs_b), known without diagonalizing."""
    import qiskit_addon_sqd
    from qiskit_addon_sqd import fermion as F
    from qiskit_addon_sqd.counts import bit_array_to_arrays
    ver = getattr(qiskit_addon_sqd, "__version__", None)
    if ver is None:
        from importlib.metadata import version
        ver = version("qiskit-addon-sqd")
    if ver != SQD_VERSION_TESTED:
        raise RuntimeError(f"sqd_batches replicates qiskit-addon-sqd {SQD_VERSION_TESTED} internals, found {ver}")
    raw_b, raw_p = bit_array_to_arrays(bit_array_from_ints(ints, norb))
    rng = np.random.default_rng(seed)
    na, nb = nelec
    md = (None, None) if max_dim is None else ((max_dim, max_dim) if isinstance(max_dim, int) else tuple(max_dim))
    cfg = F._LoopConfig(raw_bitstrings=raw_b, raw_probs=raw_p, n_alpha=na, n_beta=nb,
                        samples_per_batch=samples_per_batch, num_batches=num_batches, norb=norb,
                        symmetrize_spin=symmetrize_spin, include_a=np.array([], dtype=int),
                        include_b=np.array([], dtype=int), max_dim_a=md[0], max_dim_b=md[1],
                        energy_tol=LIN_SQD["energy_tol"], occupancies_tol=LIN_SQD["occupancies_tol"],
                        carryover_threshold=LIN_SQD["carryover_threshold"], rng=rng)
    empty = np.array([], dtype=np.int64)
    return F._prepare_ci_strings(cfg, None, empty, empty)


def solve_batch(one_body, two_body, constant, ci_strs, norb, nelec, **kw):
    """qiskit_addon_sqd.fermion.solve_sci on one batch (pyscf selected CI in the fixed space strs_a x strs_b).
    Returns (E incl. constant, info)."""
    from qiskit_addon_sqd.fermion import solve_sci
    t0 = time.time()
    res = solve_sci(ci_strs, np.asarray(one_body), np.asarray(two_body), norb=norb, nelec=tuple(nelec), **kw)
    sa, sb = ci_strs
    return float(res.energy + constant), {"dim_a": int(len(sa)), "dim_b": int(len(sb)),
                                          "dim": int(len(sa)) * int(len(sb)), "t": time.time() - t0,
                                          "spin_sq": float(res.sci_state.spin_square())}


def run_sqd(one_body, two_body, constant, ints, norb, nelec, seed=0, **overrides):
    """The library call of Lin et al. (diagonalize_fermionic_hamiltonian with LIN_SQD and solve_sci_batch) on integer
    bitstrings in the basis of (one_body, two_body).  Returns (E incl. constant, per-batch [(E, dim)])."""
    from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci_batch
    st = dict(LIN_SQD)
    st.update(overrides)
    hist = []

    def cb(results):
        hist.append([(float(r.energy + constant), int(np.prod(r.sci_state.amplitudes.shape))) for r in results])

    res = diagonalize_fermionic_hamiltonian(np.asarray(one_body), np.asarray(two_body), bit_array_from_ints(ints, norb),
                                            norb=norb, nelec=tuple(nelec), sci_solver=solve_sci_batch,
                                            seed=np.random.default_rng(seed), callback=cb, **st)
    return float(res.energy + constant), hist


# ------------------------------------------------------------------------------------------ sample statistics

def unique_counts(ints: np.ndarray, norb: int) -> dict:
    """Unique full configurations, alpha strings, beta strings and their union (the symmetrize_spin pool)."""
    a, b = ints_to_strings(ints, norb)
    ua, ub = np.unique(a), np.unique(b)
    return {"n": int(len(ints)), "unique": int(len(np.unique(ints))), "unique_alpha": int(len(ua)),
            "unique_beta": int(len(ub)), "unique_union": int(len(np.union1d(ua, ub)))}
