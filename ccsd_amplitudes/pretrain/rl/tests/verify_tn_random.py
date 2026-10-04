#!/usr/bin/env python3
"""[Independent verifier] LUCJEnergyTN vs ffsim exact_energy on small random systems (no truncation).

Covers: nocc != nvirt, nocc = 1, nvirt = 1, n_reps 1/2/3, t1 None / given, CPU (complex128) and CUDA (complex128 and
complex64) engines, basis 'er' / explicit random split basis / identity basis, non-unitary input U (polar factor),
real-valued U, large Z (phase wrap), Z = 0, and a non-symmetric Z (convention difference check only).
Usage: CUDA_VISIBLE_DEVICES=4 python3 pretrain/rl/tests/verify_tn_random.py
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ[_v] = "2"
import sys  # noqa: E402
import tempfile  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import ffsim  # noqa: E402

from pretrain.rl import tn_energy as T  # noqa: E402
from pretrain.rl.energy import exact_energy, make_ucj_op  # noqa: E402

SCR = os.environ.get("TN_VERIFY_SCRATCH", tempfile.gettempdir())


def random_unitary(n, rng):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(X)
    return Q * (np.diagonal(R) / abs(np.diagonal(R)))


def polar(U):
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def run(n, nelec, n_reps, seed, *, t1_scale=0.15, zscale=0.7, device="cpu", dtype=None, basis="er",
        nonunitary=0.0, real_u=False, z_override=None, nonsym=False, chem_ham=None, label=""):
    rng = np.random.default_rng(seed)
    if chem_ham is None:
        ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed + 1000, dtype=float)
    else:
        ham = chem_ham
    U = np.stack([random_unitary(n, rng) for _ in range(n_reps)])
    if real_u:
        U = np.stack([np.linalg.qr(rng.normal(size=(n, n)))[0] for _ in range(n_reps)]).astype(np.float64)
    Z = rng.normal(scale=zscale, size=(n_reps, n, n))
    if not nonsym:
        Z = 0.5 * (Z + Z.transpose(0, 2, 1))
    if z_override is not None:
        Z = z_override(Z)
    if nonunitary:
        U = U + nonunitary * (rng.normal(size=U.shape) + 1j * rng.normal(size=U.shape))
    t1 = rng.normal(scale=t1_scale, size=(nelec[0], n - nelec[0])) if t1_scale else None
    Uref = np.stack([polar(u) for u in U]).astype(complex)
    E_ref = exact_energy(ham, n, nelec, make_ucj_op(Z, Uref, "square", t1=t1))
    if basis == "random":
        S, occ = T.random_split_basis(n, nelec[0], rng)
        bas = (S, occ)
    elif basis == "identity":
        bas = (np.eye(n), np.arange(n) < nelec[0])
    else:
        bas = basis
    kw = {}
    if dtype is not None:
        kw["dtype"] = dtype
    ev = T.LUCJEnergyTN(ham.one_body_tensor.real, ham.two_body_tensor.real, ham.constant, n, nelec,
                        max_bond=10 ** 6, device=device, basis=bas, block2_threads=1, stack_mem=1 << 28,
                        scratch=tempfile.mkdtemp(prefix="vtn_", dir=SCR), **kw)
    E, info = ev.energy(U, Z, t1)
    # also: no-cutoff variant (cutoff=0) for the strict comparison
    ev0 = T.LUCJEnergyTN(ham.one_body_tensor.real, ham.two_body_tensor.real, ham.constant, n, nelec,
                         max_bond=10 ** 6, device=device, basis=bas, block2_threads=1, stack_mem=1 << 28,
                         cutoff=0.0, scratch=tempfile.mkdtemp(prefix="vtn_", dir=SCR), **kw)
    E0, _ = ev0.energy(U, Z, t1)
    d, d0 = E - E_ref, E0 - E_ref
    print(f"{label:28s} n={n} nelec={nelec} reps={n_reps} dev={device} dt={str(dtype).replace('torch.', '')} "
          f"basis={basis if isinstance(basis, str) else 'tuple'}: E_ref {E_ref:.10f}  dE(default) {d:+.2e}  "
          f"dE(cutoff0) {d0:+.2e}  maxbond {info['max_bond']} imag {info['imag']:+.1e}", flush=True)
    return d, d0


def chem_hamiltonian(atom, basis="sto-3g", ncore=0):
    from pyscf import gto, scf
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    md = ffsim.MolecularData.from_scf(mf, active_space=list(range(ncore, mol.nao)))
    nelec = md.nelec
    return md.hamiltonian, md.norb, nelec


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = []
    if which in ("all", "cpu"):
        res.append(("cpu nocc>nvirt", run(7, (4, 4), 2, 1, label="nocc>nvirt")))
        res.append(("cpu nocc<nvirt", run(7, (2, 2), 2, 2, label="nocc<nvirt")))
        res.append(("cpu nocc=1", run(6, (1, 1), 2, 3, label="nocc=1")))
        res.append(("cpu nvirt=1", run(6, (5, 5), 2, 4, label="nvirt=1")))
        res.append(("cpu reps=1", run(6, (3, 3), 1, 5, label="reps=1")))
        res.append(("cpu reps=3", run(6, (2, 2), 3, 6, label="reps=3")))
        res.append(("cpu t1=None", run(7, (3, 3), 2, 7, t1_scale=0.0, label="t1=None")))
        res.append(("cpu random basis", run(7, (4, 4), 2, 8, basis="random", label="random split basis")))
        res.append(("cpu identity basis", run(7, (4, 4), 2, 9, basis="identity", label="identity basis")))
        res.append(("cpu non-unitary U", run(6, (3, 3), 2, 10, nonunitary=0.05, label="non-unitary U (polar)")))
        res.append(("cpu real U", run(6, (3, 3), 2, 11, real_u=True, label="real U")))
        res.append(("cpu big Z", run(6, (3, 3), 2, 12, zscale=4.0, label="big Z (|z|>pi)")))
        res.append(("cpu Z=0", run(6, (3, 3), 2, 13, z_override=lambda Z: 0 * Z, label="Z=0")))
        res.append(("cpu big t1", run(7, (3, 3), 2, 14, t1_scale=0.8, label="big t1")))
        res.append(("cpu n=8 (5,5)", run(8, (5, 5), 2, 15, label="n8 (5,5)")))
    if which in ("all", "chem"):
        # chemistry Hamiltonian: N2 STO-3G, 2 frozen core -> 8 orbitals, (5,5)
        ham, n, ne = chem_hamiltonian("N 0 0 0; N 0 0 1.15", ncore=2)
        res.append(("cpu N2 chem", run(n, ne, 2, 16, chem_ham=ham, label="N2 sto-3g fc")))
        ham, n, ne = chem_hamiltonian("O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59", ncore=1)
        res.append(("cpu H2O chem", run(n, ne, 2, 17, chem_ham=ham, label="H2O sto-3g fc")))
        # non-symmetric Z: ffsim uses upper triangle for aa; report only
        try:
            run(6, (3, 3), 2, 18, nonsym=True, label="NONSYM Z (info only)")
        except Exception as e:  # noqa: BLE001
            print("nonsym Z raised", type(e).__name__, e)
    if which in ("all", "gpu") and torch.cuda.is_available():
        res.append(("gpu c128 nocc>nvirt", run(7, (4, 4), 2, 21, device="cuda", dtype=torch.complex128,
                                               label="gpu c128 nocc>nvirt")))
        res.append(("gpu c128 nvirt=1", run(6, (5, 5), 2, 22, device="cuda", dtype=torch.complex128,
                                            label="gpu c128 nvirt=1")))
        res.append(("gpu c128 random basis", run(7, (3, 3), 3, 23, device="cuda", dtype=torch.complex128,
                                                 basis="random", label="gpu c128 random basis")))
        ham, n, ne = chem_hamiltonian("N 0 0 0; N 0 0 1.15", ncore=2)
        res.append(("gpu c64 N2", run(n, ne, 2, 24, device="cuda", chem_ham=ham, label="gpu c64 N2")))
        res.append(("gpu c64 nocc>nvirt", run(7, (4, 4), 2, 21, device="cuda", label="gpu c64 nocc>nvirt")))
    bad = [(k, d) for k, d in res if abs(d[1]) > 1e-8 and "c64" not in k]
    print("worst |dE| (cutoff 0, c128):", max((abs(d[1]) for k, d in res if "c64" not in k), default=0))
    print("worst |dE| (default cutoff, c128):", max((abs(d[0]) for k, d in res if "c64" not in k), default=0))
    print("FAILURES:" if bad else "no failures (1e-8)", bad)
