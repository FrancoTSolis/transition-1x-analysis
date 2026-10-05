#!/usr/bin/env python3
"""Tests of pretrain/rl/tn_sampler.py (MPS direct sampling for TN-sampled SQD) on small random systems.

  1. amplitudes: MPS amplitudes (dm engine / torch, and the numpy zip-up engine) == ffsim state in the sampled basis
     c = O(Fr)^dag psi = ffsim.apply_orbital_rotation(psi, Fr^T), up to ffsim_sign and one global phase;
  2. basis: <c|H^chi|c> with sampled_basis_integrals == ffsim MO-basis energy; mps_energy == ev.energy;
  3. sampler: empirical frequencies of MPSSampler.sample == |c|^2 (chi-square test), logp == log|amp|^2, also for a
     truncated MPS (then against the truncated MPS's own amplitudes);
  4. integer bitstrings: occ_to_ints follows ffsim BitstringType.INT (frequencies vs ffsim.sample_state_vector);
  5. exact references (CUDA): LUCJStateGPU(basis=Fr) CI matrix == c, energies equal; sample_ci_matrix distribution;
  6. SQD: min_b solve_batch(sqd_batches(...)[b]) == qiskit_addon_sqd diagonalize_fermionic_hamiltonian (same seed,
     max_iterations 1); full string space -> FCI energy, identical in the MO and the sampled basis.
Usage: CUDA_VISIBLE_DEVICES=<gpu> python3 pretrain/rl/tests/test_tn_sampler.py   (CPU-only parts run without a GPU)
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

try:                                   # small test systems -> numba "grid size" warnings from gpu_energy
    from numba.core.errors import NumbaPerformanceWarning
    warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)
except Exception:  # noqa: BLE001
    pass

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_energy as T  # noqa: E402
from pretrain.rl import tn_sampler as S  # noqa: E402

SCR = tempfile.mkdtemp(prefix="tn_sampler_test_", dir=os.environ.get("TN_SCRATCH", None))


def random_unitary(n, rng):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(X)
    return Q * (np.diagonal(R) / abs(np.diagonal(R)))


def system(n, k, seed, zscale=0.5, t1scale=0.3):
    import ffsim
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U = np.stack([random_unitary(n, rng) for _ in range(2)])
    Z = rng.normal(scale=zscale, size=(2, n, n))
    Z = Z + Z.transpose(0, 2, 1)
    t1 = rng.normal(scale=t1scale, size=(k, n - k))
    Sb, occ = T.random_split_basis(n, k, rng)
    return ham, U, Z, t1, Sb, occ


def ffsim_state(ham, n, k, U, Z, t1):
    import ffsim
    from pretrain.rl.energy import make_ucj_op
    return ffsim.apply_unitary(ffsim.hartree_fock_state(n, (k, k)), make_ucj_op(Z, U, "square", t1=t1),
                               norb=n, nelec=(k, k))


def evaluator(ham, n, k, Sb, occ, device="cpu", dtype=torch.complex128, chi=10 ** 6, cutoff=0.0, method="dm"):
    return T.LUCJEnergySplitTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, (k, k), Sb, occ,
                               max_bond=chi, cutoff=cutoff, device=device, dtype=dtype, block2_threads=1,
                               scratch=tempfile.mkdtemp(dir=SCR), stack_mem=1 << 28, method=method)


def all_occ(n, k):
    strs = S.ci_strings(n, k)
    a = np.repeat(strs, len(strs))
    b = np.tile(strs, len(strs))
    return S.strings_to_occ(a, b, n), a, b             # row order = ffsim (alpha address major)


def chi2_pvalue(counts, probs, n):
    """Pearson chi-square goodness of fit (bins with expected >= 5, the rest pooled)."""
    from scipy.stats import chi2
    e = probs * n
    big = e >= 5
    obs = np.r_[counts[big], counts[~big].sum()]
    exp_ = np.r_[e[big], e[~big].sum()]
    keep = exp_ > 0
    stat = float(((obs[keep] - exp_[keep]) ** 2 / exp_[keep]).sum())
    return float(chi2.sf(stat, keep.sum() - 1)), stat, int(keep.sum() - 1)


def test_amplitudes(device="cpu"):
    import ffsim
    for n, k, seed, method in ((5, 2, 1, "dm"), (6, 3, 2, "dm"), (7, 3, 4, "dm"), (7, 2, 5, "dm"), (6, 3, 7, "zipup")):
        if method == "zipup" and device != "cpu":
            continue
        ham, U, Z, t1, Sb, occ = system(n, k, seed)
        psi = ffsim_state(ham, n, k, U, Z, t1)
        ev = evaluator(ham, n, k, Sb, occ, device=device, method=method)
        mps, Fr, info = S.build_mps(ev, U, Z, t1)
        assert np.allclose(Fr, S.sampled_basis(ev, t1))
        c = ffsim.apply_orbital_rotation(psi, Fr.T, norb=n, nelec=(k, k))
        smp = S.MPSSampler(mps)
        O, _, _ = all_occ(n, k)
        amp = smp.amplitudes(O) * S.ffsim_sign(O)
        j = np.argmax(np.abs(c))
        ph = c[j] / amp[j]
        err = np.abs(c - ph * amp).max()
        print(f"  amplitudes {method} {device} n={n} k={k}: max|c - phase*sign*amp| {err:.1e} |phase| {abs(ph):.12f} "
              f"isometry err {smp.isometry_err:.1e}")
        assert err < 1e-10 and abs(abs(ph) - 1) < 1e-10
        # wrong-sector configuration -> zero amplitude
        bad = O[:1].copy()
        bad[0, np.argmax(bad[0] == 0)] = 1
        assert abs(smp.amplitudes(bad)[0]) == 0.0


def test_basis_energy():
    import ffsim
    for n, k, seed in ((6, 3, 11), (7, 2, 12)):
        ham, U, Z, t1, Sb, occ = system(n, k, seed)
        psi = ffsim_state(ham, n, k, U, Z, t1)
        E0 = float(np.vdot(psi, ffsim.linear_operator(ham, n, (k, k)) @ psi).real)
        ev = evaluator(ham, n, k, Sb, occ)
        mps, Fr, _ = S.build_mps(ev, U, Z, t1)
        h, g = S.sampled_basis_integrals(ham.one_body_tensor, ham.two_body_tensor, Fr)
        assert h.dtype == np.float64 and g.dtype == np.float64
        c = ffsim.apply_orbital_rotation(psi, Fr.T, norb=n, nelec=(k, k))
        hamF = ffsim.MolecularHamiltonian(h, g, ham.constant)
        E1 = float(np.vdot(c, ffsim.linear_operator(hamF, n, (k, k)) @ c).real)
        Em, _ = S.mps_energy(ev, mps, Fr)
        Ee, _ = ev.energy(U, Z, t1)
        print(f"  basis n={n}: E_MO {E0:.10f}  E_chi(ffsim) {E1 - E0:+.1e}  mps_energy {Em - E0:+.1e}  "
              f"ev.energy {Ee - E0:+.1e}")
        assert abs(E1 - E0) < 1e-9 and abs(Em - E0) < 1e-8 and abs(Ee - Em) < 1e-10


def test_sampler(device="cpu"):
    for n, k, seed, chi in ((6, 3, 21, 10 ** 6), (7, 3, 22, 10 ** 6), (8, 4, 23, 12)):
        ham, U, Z, t1, Sb, occ = system(n, k, seed, zscale=1.0)
        ev = evaluator(ham, n, k, Sb, occ, device=device, chi=chi, cutoff=0.0 if chi > 1000 else 1e-12)
        mps, Fr, info = S.build_mps(ev, U, Z, t1)
        smp = S.MPSSampler(mps)
        O, _, _ = all_occ(n, k)
        p = smp.probabilities(O)
        N = 200_000
        occ_s, lp = smp.sample(N, seed=seed, batch=65_536)
        ints_all = S.occ_to_ints(O)
        ints_s = S.occ_to_ints(occ_s)
        assert np.isin(ints_s, ints_all).all(), "sample outside the particle-number sector"
        order = np.argsort(ints_all)
        cnt = np.bincount(np.searchsorted(ints_all[order], ints_s), minlength=len(ints_all))
        cnt_ffsim_order = np.empty_like(cnt)
        cnt_ffsim_order[order] = cnt
        pv, stat, dof = chi2_pvalue(cnt_ffsim_order, p, N)
        lp_ref = np.log(smp.probabilities(occ_s[:2000]))
        dlp = np.abs(lp[:2000] - lp_ref).max()
        tvd = 0.5 * np.abs(cnt_ffsim_order / N - p).sum()
        print(f"  sampler {device} n={n} k={k} chi={chi} (max bond {info['max_bond']}, disc {info['discarded_sum']:.1e}): "
              f"sum p {p.sum():.12f}  chi2 p-value {pv:.3f} (stat {stat:.0f}, dof {dof})  TVD {tvd:.4f}  "
              f"max|logp - log|amp|^2| {dlp:.1e}")
        assert abs(p.sum() - 1) < 1e-10 and pv > 1e-4 and dlp < 1e-8


def test_int_convention():
    """occ_to_ints == ffsim BitstringType.INT: identical distributions of ffsim.sample_state_vector(c) and the MPS."""
    import ffsim
    n, k, seed = 6, 2, 31
    ham, U, Z, t1, Sb, occ = system(n, k, seed, zscale=1.0)
    psi = ffsim_state(ham, n, k, U, Z, t1)
    ev = evaluator(ham, n, k, Sb, occ)
    mps, Fr, _ = S.build_mps(ev, U, Z, t1)
    c = ffsim.apply_orbital_rotation(psi, Fr.T, norb=n, nelec=(k, k))
    N = 100_000
    f_ints = np.asarray(ffsim.sample_state_vector(c, norb=n, nelec=(k, k), shots=N, seed=1,
                                                  bitstring_type=ffsim.BitstringType.INT), dtype=np.int64)
    m_ints = S.occ_to_ints(S.MPSSampler(mps).sample(N, seed=2)[0])
    # exact probabilities keyed by the ffsim-convention integer (alpha = low bits)
    strs = S.ci_strings(n, k)
    keys = (strs[:, None] | (strs[None, :] << n)).reshape(-1)
    p = np.abs(c) ** 2
    pv = []
    for ints in (f_ints, m_ints):
        cnt = np.array([np.count_nonzero(ints == kk) for kk in keys])
        pv.append(chi2_pvalue(cnt, p, N)[0])
    # HF in the sampled basis: alpha = beta = sites of the occupied split orbitals
    hf = S.occ_to_ints(np.array([[3 if o else 0 for o in occ]], dtype=np.uint8))[0]
    hf_ffsim = int(sum(1 << j for j in range(n) if occ[j])) * (1 + (1 << n))
    print(f"  int convention: chi2 p-values ffsim {pv[0]:.3f}  MPS {pv[1]:.3f}; HF int {hf} == {hf_ffsim}")
    assert min(pv) > 1e-4 and hf == hf_ffsim


def test_exact_gpu():
    if not torch.cuda.is_available():
        print("  (no CUDA: LUCJStateGPU / sample_ci_matrix skipped)")
        return
    import ffsim
    for n, k, seed in ((6, 3, 41), (8, 3, 42), (9, 4, 43)):
        ham, U, Z, t1, Sb, occ = system(n, k, seed)
        psi = ffsim_state(ham, n, k, U, Z, t1)
        E0 = float(np.vdot(psi, ffsim.linear_operator(ham, n, (k, k)) @ psi).real)
        Fr = S.sampled_basis(type("ev", (), {"S": Sb})(), t1)
        h, g = S.sampled_basis_integrals(ham.one_body_tensor, ham.two_body_tensor, Fr)
        c = ffsim.apply_orbital_rotation(psi, Fr.T, norb=n, nelec=(k, k))
        for dt, tol_v, tol_e in ((torch.complex128, 1e-10, 1e-9), (torch.complex64, 2e-5, 2e-4)):
            eng = S.LUCJStateGPU(h, g, ham.constant, n, (k, k), basis=Fr, dtype=dt, max_mem_gb=1.0)
            E = eng.energy(U, Z, t1)
            v = eng.state(U, Z, t1).reshape(-1).cpu().numpy().astype(np.complex128)
            j = np.argmax(np.abs(c))
            ph = c[j] / v[j]
            err = np.abs(c - ph * v).max()
            print(f"  LUCJStateGPU {str(dt)[6:]} n={n} k={k}: max|c - phase v| {err:.1e}  E_chi - E_MO {E - E0:+.1e}")
            assert err < tol_v and abs(E - E0) < tol_e * max(1.0, abs(E0))
            eng.release()
        # exact sampler on the GPU
        eng = S.LUCJStateGPU(h, g, ham.constant, n, (k, k), basis=Fr, dtype=torch.complex128, max_mem_gb=1.0)
        P = eng.state(U, Z, t1).clone()
        eng.release()
        N = 200_000
        ia, ib = S.sample_ci_matrix(P, N, seed=5)
        dim = P.shape[1]
        cnt = np.bincount(ia * dim + ib, minlength=dim * dim)
        pv = chi2_pvalue(cnt, np.abs(c) ** 2, N)[0]
        print(f"  sample_ci_matrix n={n}: chi2 p-value {pv:.3f}")
        assert pv > 1e-4


def test_sqd():
    """sqd_batches + solve_batch reproduce the library call; full string space = FCI in both bases."""
    import ffsim
    from pyscf import fci
    n, k, seed = 8, 3, 51
    ham, U, Z, t1, Sb, occ = system(n, k, seed, zscale=1.0)
    psi = ffsim_state(ham, n, k, U, Z, t1)
    ev = evaluator(ham, n, k, Sb, occ)
    mps, Fr, _ = S.build_mps(ev, U, Z, t1)
    h, g = S.sampled_basis_integrals(ham.one_body_tensor, ham.two_body_tensor, Fr)
    ints = S.occ_to_ints(S.MPSSampler(mps).sample(3000, seed=3)[0])
    kw = dict(samples_per_batch=60, num_batches=3, max_dim=24)
    batches = S.sqd_batches(ints, n, (k, k), seed=7, **kw)
    Eb = [S.solve_batch(h, g, ham.constant, b, n, (k, k))[0] for b in batches]
    E_lib, hist = S.run_sqd(h, g, ham.constant, ints, n, (k, k), seed=7, **kw)
    dims = [d for _, d in hist[0]]
    print(f"  sqd: batches {[f'{e:.8f}' for e in Eb]} dims {[len(a) * len(b) for a, b in batches]} | library "
          f"{E_lib:.8f} dims {dims}")
    assert abs(min(Eb) - E_lib) < 1e-9 and dims == [len(a) * len(b) for a, b in batches]
    # all strings -> FCI, the same in the MO basis and in the sampled basis
    allstr = S.ci_strings(n, k)
    # (random Hamiltonian, |E| ~ 1e3 Ha: tight Davidson tolerances in both solvers)
    e_mo = fci.direct_spin1.kernel(ham.one_body_tensor, ham.two_body_tensor, n, (k, k), tol=1e-14,
                                   max_cycle=500)[0] + ham.constant
    e_chi = S.solve_batch(h, g, ham.constant, (allstr, allstr), n, (k, k), tol=1e-14, max_cycle=500)[0]
    print(f"  sqd full space: FCI(MO) {e_mo:.10f}  QSCI(all strings, chi basis) {e_chi:.10f}  diff {e_chi - e_mo:+.1e}")
    assert abs(e_chi - e_mo) < 1e-9 * abs(e_mo)


if __name__ == "__main__":
    devs = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
    for dev in devs:
        test_amplitudes(dev)
        test_sampler(dev)
    test_basis_energy()
    test_int_convention()
    test_exact_gpu()
    test_sqd()
    print("all ok")
