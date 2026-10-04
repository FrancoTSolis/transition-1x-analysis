#!/usr/bin/env python3
"""Correctness tests of pretrain/rl/tn_energy.py on tiny random systems (dense references, CPU, complex128).

  1. Givens sequence reconstructs the orbital rotation;
  2. LUCJ MPS (no truncation) energy == ffsim exact energy (random Hamiltonian, random U, Z, t1);
  3. the MPS state itself == ffsim state up to the fermionic sign/permutation map (via energies of several H).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pretrain.rl import tn_energy as T  # noqa: E402


def random_unitary(n, rng):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(X)
    return Q * (np.diagonal(R) / abs(np.diagonal(R)))


def test_givens(n=6, seed=0):
    rng = np.random.default_rng(seed)
    R = random_unitary(n, rng)
    seq, ph = T.givens_sequence(R)
    M = np.eye(n, dtype=complex)
    for lo, g in seq:
        E = np.eye(n, dtype=complex)
        E[lo:lo + 2, lo:lo + 2] = g
        M = E @ M
    M = np.diag(ph) @ M
    err = np.abs(M - R).max()
    print(f"givens reconstruction err {err:.2e}")
    assert err < 1e-10


def random_params(n, rng, n_reps=2, zscale=0.5):
    U = np.stack([random_unitary(n, rng) for _ in range(n_reps)])
    Z = rng.normal(scale=zscale, size=(n_reps, n, n))
    Z = Z + Z.transpose(0, 2, 1)
    return U, Z


def test_energy(n=5, nelec=(2, 2), seed=1, with_t1=True):
    import ffsim
    from pretrain.rl.energy import exact_energy, make_ucj_op
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U, Z = random_params(n, rng)
    t1 = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0])) if with_t1 else None
    E_ref = exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=t1))
    rots, dcs, F = T.lucj_layers(U, Z, t1)
    for slater in ((False, True) if nelec[0] == nelec[1] else (False,)):
        mps = T.SymMPS(n, nelec, device="cpu")
        for k, (R, (znn, zon)) in enumerate(zip(rots, dcs)):
            if k == 0 and slater:
                mps.apply_gate_sequence(T.slater_sequence(R, nelec[0]), max_bond=10 ** 6)
            else:
                mps.apply_orbital_rotation(R, max_bond=10 ** 6)
            mps.apply_diag_coulomb(znn, zon, max_bond=10 ** 6)
        v = mps.to_dense()
        h, g = T.rotate_hamiltonian(ham.one_body_tensor, ham.two_body_tensor, F)
        H = T.dense_fock_hamiltonian(h, g, ham.constant, n)
        E = float(np.vdot(v, H @ v).real / np.vdot(v, v).real)
        print(f"n={n} nelec={nelec} slater={slater}: E_mps {E:.10f}  E_ffsim {E_ref:.10f}  diff {E - E_ref:.2e} "
              f"| norm {np.vdot(v, v).real:.12f} bonds {mps.bond_dims()} disc {mps.discarded:.1e}")
        assert abs(E - E_ref) < 1e-8


def test_block2(n=5, nelec=(2, 2), seed=1, sign_ab=1.0):
    import ffsim
    import tempfile
    from pretrain.rl.energy import exact_energy, make_ucj_op
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U, Z = random_params(n, rng)
    t1 = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
    E_ref = exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=t1))
    rots, dcs, F = T.lucj_layers(U, Z, t1)
    mps = T.SymMPS(n, nelec, device="cpu")
    for R, (znn, zon) in zip(rots, dcs):
        mps.apply_orbital_rotation(R, max_bond=10 ** 6)
        mps.apply_diag_coulomb(znn, zon, max_bond=10 ** 6)
    h, g = T.rotate_hamiltonian(ham.one_body_tensor, ham.two_body_tensor, F)
    tmp = tempfile.mkdtemp(prefix="b2t_", dir=os.environ.get("TN_SCRATCH", "/tmp"))
    B = T.Block2Energy(n, nelec, tmp, n_threads=2, stack_mem=1 << 28, sign_ab=sign_ab)
    bm = B.to_block2(mps)
    mpo = B.mpo(h, g, ham.constant)
    E = B.expectation(bm, mpo) / B.norm2(bm)
    print(f"block2 n={n} nelec={nelec} sign_ab={sign_ab}: E {E.real:.10f} (imag {E.imag:.1e}) ffsim {E_ref:.10f} "
          f"diff {E.real - E_ref:.2e}")
    return E.real - E_ref


def test_fishman_white(n=8, nelec=(3, 3), seed=5, window=8):
    import ffsim
    from pretrain.rl.energy import exact_energy, make_ucj_op
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U, Z = random_params(n, rng)
    E_ref = exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=None))
    rots, dcs, F = T.lucj_layers(U, Z, None)
    occ, seq, infid = T.fishman_white_sequence(rots[0][:, :nelec[0]], window=window)
    mps = T.SymMPS(n, nelec, device="cpu", occ=[3 * int(o) for o in occ])
    mps.apply_gate_sequence(seq, max_bond=10 ** 6)
    mps.apply_diag_coulomb(*dcs[0], max_bond=10 ** 6)
    mps.apply_orbital_rotation(rots[1], max_bond=10 ** 6)
    mps.apply_diag_coulomb(*dcs[1], max_bond=10 ** 6)
    v = mps.to_dense()
    h, g = T.rotate_hamiltonian(ham.one_body_tensor, ham.two_body_tensor, F)
    H = T.dense_fock_hamiltonian(h, g, ham.constant, n)
    E = float(np.vdot(v, H @ v).real / np.vdot(v, v).real)
    print(f"FW n={n} window={window}: infid {infid:.2e} occ {occ.tolist()} E {E:.10f} ref {E_ref:.10f} diff {E - E_ref:.2e}")
    if window >= n:
        assert abs(E - E_ref) < 1e-8


def test_split(n=6, nelec=(3, 3), seed=11, max_bond=10 ** 6):
    import ffsim
    import tempfile
    from pretrain.rl.energy import exact_energy, make_ucj_op
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U, Z = random_params(n, rng)
    t1 = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
    E_ref = exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=t1))
    S, occ = T.random_split_basis(n, nelec[0], rng)
    tmp = tempfile.mkdtemp(prefix="b2s_", dir=os.environ.get("TN_SCRATCH", "/tmp"))
    ev = T.LUCJEnergySplitTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, nelec, S, occ,
                             max_bond=max_bond, cutoff=0.0, device="cpu", dtype=__import__("torch").complex128,
                             block2_threads=2, scratch=tmp, stack_mem=1 << 28)
    E, info = ev.energy(U, Z, t1)
    print(f"split n={n} nelec={nelec} chi={max_bond}: E {E:.10f} ffsim {E_ref:.10f} diff {E - E_ref:.2e} "
          f"maxbond {info['max_bond']} disc {info['discarded_sum']:.1e}")
    if max_bond >= 10 ** 5:
        assert abs(E - E_ref) < 1e-8


def test_public_api(n=7, nelec=(3, 3), seed=21):
    """LUCJEnergyTN with the integral-only ER split basis, no truncation, CPU engine: equals ffsim exactly."""
    import ffsim
    import tempfile
    from pretrain.rl.energy import exact_energy, make_ucj_op
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U, Z = random_params(n, rng, zscale=0.3)
    t1 = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
    E_ref = exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=t1))
    ev = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, nelec, max_bond=10 ** 6,
                        device="cpu", basis="er", cutoff=0.0, block2_threads=1, stack_mem=1 << 28, zip_margin=1.0,
                        scratch=tempfile.mkdtemp(prefix="b2p_", dir=os.environ.get("TN_SCRATCH", "/tmp")))
    E, info = ev.energy(U, Z, t1)
    print(f"public API (ER basis, cpu) n={n}: E {E:.10f} ffsim {E_ref:.10f} diff {E - E_ref:.2e} "
          f"maxbond {info['max_bond']}")
    # the zip-up keeps a 1e-15 relative noise cutoff (NpSymMPS.zip_cutoff); the LUCJ state is not an eigenstate of H,
    # so the energy is first order in the dropped amplitude: O(1e-8) on these random large-integral Hamiltonians
    assert abs(E - E_ref) < 1e-6


if __name__ == "__main__":
    import sys as _s
    if len(_s.argv) > 1 and _s.argv[1] == "block2":
        for sg in (1.0, -1.0):
            test_block2(4, (2, 2), 3, sg)
            test_block2(5, (3, 2), 2, sg)
        raise SystemExit
    test_givens()
    test_energy(4, (2, 2), seed=3)
    test_energy(5, (2, 2), seed=1)
    test_energy(5, (3, 2), seed=2)
    test_energy(6, (3, 3), seed=4)
    test_split(6, (3, 3), 11)
    test_split(7, (3, 3), 12)
    test_public_api(6, (3, 3), 21)
    print("all ok")
