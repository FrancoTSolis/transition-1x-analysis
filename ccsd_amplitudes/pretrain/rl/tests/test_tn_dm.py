#!/usr/bin/env python3
"""Correctness tests of the density-matrix (method="dm") TN energy against ffsim, plus engine invariants.

  1. no truncation: DM energy == ffsim exact_energy(make_ucj_op(Z, U, "square", t1)) on random Hamiltonians, random
     split bases, n_reps 1-3, closed shells with nocc <, =, > nvirt, t1 None / large, |Z| > pi, Z = 0
     (complex128, torch engine on the CPU; and on CUDA when available);
  2. real molecules (N2, H2O, STO-3G, frozen core): complex128 and complex64 vs ffsim;
  3. block-QR canonical moves leave the state unchanged and produce isometries;
  4. public API (LUCJEnergyTN, ER basis, both methods) without truncation; info / settings keys;
  5. with truncation: discarded weights are reported (dm: discarded_sum; zipup: also discarded_zip_sum).
Usage: python3 pretrain/rl/tests/test_tn_dm.py      (TN_MODULE=pretrain.rl.tn_energy by default)
"""
from __future__ import annotations

import importlib
import os
import sys
import tempfile
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
T = importlib.import_module(os.environ.get("TN_MODULE", "pretrain.rl.tn_energy"))
SCR = tempfile.mkdtemp(prefix="tn_dm_test_", dir=os.environ.get("TN_SCRATCH", None))


def random_unitary(n, rng):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(X)
    return Q * (np.diagonal(R) / abs(np.diagonal(R)))


def ffsim_energy(ham, n, nelec, U, Z, t1):
    from pretrain.rl.energy import exact_energy, make_ucj_op
    return exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=t1))


def split_ev(ham, n, nelec, S, occ, device, dtype, cutoff, chi=10 ** 6, method="dm"):
    return T.LUCJEnergySplitTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, nelec, S, occ,
                               max_bond=chi, cutoff=cutoff, device=device, dtype=dtype, block2_threads=1,
                               scratch=tempfile.mkdtemp(dir=SCR), stack_mem=1 << 28, method=method)


def test_random(device="cpu"):
    import ffsim
    worst = 0.0
    cases = [(4, (2, 2), 3, 2, 0.0), (5, (2, 2), 1, 2, 0.0), (6, (3, 3), 4, 2, 0.0), (6, (1, 1), 8, 2, 0.0),
             (6, (5, 5), 9, 2, 0.0), (6, (3, 3), 10, 1, 0.0), (6, (2, 2), 11, 3, 0.0),
             (7, (3, 3), 12, 2, 1e-14), (7, (4, 4), 5, 2, 1e-14), (7, (2, 2), 6, 2, 1e-14)]
    for n, nelec, seed, n_reps, cutoff in cases:
        rng = np.random.default_rng(seed)
        ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
        U = np.stack([random_unitary(n, rng) for _ in range(n_reps)])
        Z = rng.normal(scale=2.0 if seed % 2 else 0.5, size=(n_reps, n, n))       # |Z| > pi for odd seeds
        Z = Z + Z.transpose(0, 2, 1)
        if seed == 10:
            Z[:] = 0.0
        t1 = None if seed % 3 == 0 else rng.normal(scale=0.5, size=(nelec[0], n - nelec[0]))
        E_ref = ffsim_energy(ham, n, nelec, U, Z, t1)
        S, occ = T.random_split_basis(n, nelec[0], rng)
        E, info = split_ev(ham, n, nelec, S, occ, device, torch.complex128, cutoff).energy(U, Z, t1)
        d = E - E_ref
        # Gram-based truncation resolves amplitudes down to ~sqrt(eps) relative; the energy is first order in the
        # amplitude, and these random Hamiltonians have |E| ~ 10-800 Ha: 1e-7 (no cutoff) / 1e-5 (cutoff 1e-14)
        tol = 1e-7 if cutoff == 0.0 else 1e-5
        print(f"  dm {device} n={n} nelec={nelec} reps={n_reps} t1={'no' if t1 is None else 'yes'} cutoff={cutoff:g}: "
              f"E {E:.10f} ffsim {E_ref:.10f} diff {d:+.2e} maxbond {info['max_bond']} n_dm {info['n_dm']}", flush=True)
        assert abs(d) < tol, d
        worst = max(worst, abs(d) if cutoff == 0.0 else 0.0)
    print(f"random systems ({device}): worst |dE| without cutoff {worst:.2e}")


def molecule(kind):
    import ffsim
    import pyscf
    if kind == "N2":
        mol = pyscf.gto.M(atom="N 0 0 0; N 0 0 1.1", basis="sto-3g", verbose=0)
    else:
        mol = pyscf.gto.M(atom="O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", basis="sto-3g", verbose=0)
    mf = pyscf.scf.RHF(mol).run()
    data = ffsim.MolecularData.from_scf(mf, active_space=range(1, mol.nao))
    return data.hamiltonian, data.norb, data.nelec


def test_molecules(device):
    for kind in (("N2", "H2O") if device != "cpu" else ("H2O",)):      # N2 (9 orbitals) is slow on a CPU core
        ham, n, nelec = molecule(kind)
        rng = np.random.default_rng(7)
        U = np.stack([random_unitary(n, rng) for _ in range(2)])
        Z = rng.normal(scale=0.5, size=(2, n, n))
        Z = Z + Z.transpose(0, 2, 1)
        t1 = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
        E_ref = ffsim_energy(ham, n, nelec, U, Z, t1)
        S, occ = T.split_basis_from_integrals(ham.one_body_tensor, ham.two_body_tensor, nelec[0])
        for dt, cut, tol in ((torch.complex128, 1e-14, 1e-8), (torch.complex64, 1e-9, 1e-4)):
            if dt == torch.complex64 and device == "cpu":
                continue
            E, info = split_ev(ham, n, nelec, S, occ, device, dt, cut).energy(U, Z, t1)
            print(f"  {kind} {device} {str(dt)[6:]}: E {E:.10f} ffsim {E_ref:.10f} diff {E - E_ref:+.2e}", flush=True)
            assert abs(E - E_ref) < tol


def test_local_windows(device):
    """Identity basis and orbital rotations that only mix neighbouring orbitals: the factor windows are 1-4 sites,
    so the centre has to be moved between factors (block-QR moves) and both sweep directions occur."""
    import ffsim
    from scipy.linalg import expm
    for n, nelec, seed in ((6, (3, 3), 31), (7, (2, 2), 32)):
        rng = np.random.default_rng(seed)
        ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
        Us = []
        for k in range(2):
            K = np.zeros((n, n), dtype=complex)
            for p in range(k % 2, n - 1, 2):                    # disjoint neighbour pairs, alternating offsets
                a = rng.normal() + 1j * rng.normal()
                K[p, p + 1], K[p + 1, p] = a, -np.conj(a)
            Us.append(expm(K))
        U = np.stack(Us)
        Z = rng.normal(scale=0.5, size=(2, n, n))
        Z = Z + Z.transpose(0, 2, 1)
        t1 = rng.normal(scale=0.2, size=(nelec[0], n - nelec[0]))
        E_ref = ffsim_energy(ham, n, nelec, U, Z, t1)
        S, occ = np.eye(n), np.arange(n) < nelec[0]
        ev = split_ev(ham, n, nelec, S, occ, device, torch.complex128, 0.0)
        E, info = ev.energy(U, Z, t1)
        print(f"  local windows {device} n={n} nelec={nelec}: E {E:.10f} ffsim {E_ref:.10f} diff {E - E_ref:+.2e} "
              f"maxbond {info['max_bond']} n_dm {info['n_dm']} t_qr {info['t_qr']:.2f}s", flush=True)
        assert abs(E - E_ref) < 1e-7


def test_moves(device):
    import ffsim
    n, nelec = 6, (3, 3)
    rng = np.random.default_rng(3)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=3, dtype=float)
    U = np.stack([random_unitary(n, rng) for _ in range(2)])
    Z = rng.normal(scale=0.5, size=(2, n, n))
    Z = Z + Z.transpose(0, 2, 1)
    S, occ = T.random_split_basis(n, nelec[0], rng)
    ev = split_ev(ham, n, nelec, S, occ, device, torch.complex128, 0.0)
    mps, _ = ev.state(U, Z, None)
    v0 = mps.to_dense()
    for target in (n - 1, 0, 3):
        mps.move_center(target)
        v = mps.to_dense()
        err = np.abs(v - v0).max()
        iso = 0.0
        for i in range(n):
            A = mps.A[i].to(torch.complex128)
            if i < mps.center:
                G = torch.einsum("asb,asc->bc", A.conj(), A)
            elif i > mps.center:
                G = torch.einsum("asb,csb->ac", A.conj(), A)
            else:
                continue
            iso = max(iso, float((G - torch.eye(G.shape[0], dtype=G.dtype, device=G.device)).abs().max()))
        print(f"  move_center({target}) {device}: max|psi - psi0| {err:.1e}, isometry error {iso:.1e}, "
              f"bonds {mps.bond_dims()}")
        assert err < 1e-12 and iso < 1e-12


def test_public_api():
    import ffsim
    n, nelec = 7, (3, 3)
    rng = np.random.default_rng(21)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=21, dtype=float)
    U = np.stack([random_unitary(n, rng) for _ in range(2)])
    Z = rng.normal(scale=0.3, size=(2, n, n))
    Z = Z + Z.transpose(0, 2, 1)
    t1 = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
    E_ref = ffsim_energy(ham, n, nelec, U, Z, t1)
    for method, kw in (("dm", {"cutoff": 0.0}), ("zipup", {"cutoff": 0.0, "zip_margin": 1.0})):
        ev = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, nelec, max_bond=10 ** 6,
                            device="cpu", basis="er", block2_threads=1, stack_mem=1 << 28, method=method,
                            scratch=tempfile.mkdtemp(dir=SCR), **kw)
        E, info = ev.energy(U, Z, t1)
        st = ev.settings()
        print(f"  public API {method} (ER basis, cpu): E {E:.10f} ffsim {E_ref:.10f} diff {E - E_ref:+.2e}  "
              f"settings {st}")
        for key in ("method", "chi", "discarded_sum", "discarded_max", "discarded_zip_sum", "discarded_total",
                    "max_bond", "t_state", "t_expect", "t_total"):
            assert key in info, key
        for key in ("method", "zip_margin", "cutoff", "dtype", "device", "engine", "basis", "mode_tol"):
            assert key in st, key
        assert abs(E - E_ref) < (1e-9 if method == "dm" else 1e-6)
    # default method: "dm" on CUDA, "zipup" on the CPU (performance; see LUCJEnergyTN)
    for dev in ["cpu"] + (["cuda"] if torch.cuda.is_available() else []):
        ev = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, nelec, max_bond=16,
                            device=dev, basis="er", stack_mem=1 << 28, scratch=tempfile.mkdtemp(dir=SCR))
        print(f"  default method on {dev}: {ev.settings()['method']}")
        assert ev.settings()["method"] == ("zipup" if dev == "cpu" else "dm")


def test_truncation_report(device):
    import ffsim
    n, nelec = 8, (4, 4)
    rng = np.random.default_rng(5)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=5, dtype=float)
    U = np.stack([random_unitary(n, rng) for _ in range(2)])
    Z = rng.normal(scale=0.5, size=(2, n, n))
    Z = Z + Z.transpose(0, 2, 1)
    S, occ = T.random_split_basis(n, nelec[0], rng)
    E_ref = ffsim_energy(ham, n, nelec, U, Z, None)
    for method in ("dm", "zipup"):
        ev = split_ev(ham, n, nelec, S, occ, device, torch.complex128, 1e-12, chi=12, method=method)
        E, info = ev.energy(U, Z, None)
        print(f"  chi 12 {method} {device}: err {E - E_ref:+.3e}  discarded_sum {info['discarded_sum']:.2e} "
              f"zip {info['discarded_zip_sum']:.2e} maxbond {info['max_bond']}")
        assert info["max_bond"] <= 12 and info["discarded_sum"] > 0
        if method == "zipup":
            assert info["discarded_zip_sum"] > 0


if __name__ == "__main__":
    devs = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
    for dev in devs:
        test_random(dev)
        test_molecules(dev)
        test_moves(dev)
        test_local_windows(dev)
        test_truncation_report(dev)
    test_public_api()
    print("all ok")
