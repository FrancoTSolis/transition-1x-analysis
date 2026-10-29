"""Shared helpers for running Avijit's UCJ-BCH2 optimizer (avijit_simulator/UCJ_optimizer.py, unmodified).

* `load_avijit()` imports the optimizer module from its file (the folder's __init__.py imports a
  `pyscf.ucj` package that does not exist, so the folder cannot be imported as a package).
* `small_mf(name)` builds the RHF objects of the test systems named in UCJ_context.md (H4, LiH, C2).
* `active_mf(h1, eri, nocc, const)` wraps an active-space Hamiltonian (our frozen-core rhf_hamiltonians/*.npz)
  in a fake-molecule RHF, so that the optimizer runs unchanged in the same orbital space as the network.
* `ffsim_energy(...)` is the independent exact energy of Avijit's ansatz.  The mapping (verified in
  `bch2_reproduce.py`): the optimizer's state is exp(i X)|HF> with X = R(G) for G = 1/2 sum g_pq n_p n_q and
  R the orbital rotation `rotation`; in ffsim this is UCJOpSpinBalanced(n_reps=1,
  diag_coulomb_mats=[[g, g]], orbital_rotations=[rotation.T]) (ffsim applies U exp(iJ) U^dag with U's columns
  the new orbitals; the optimizer's `_Operator.rotate(R)` uses R's rows).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]          # ccsd_amplitudes/
AVIJIT = ROOT / "avijit_simulator" / "UCJ_optimizer.py"


class _Ao2moShim:
    """Stands in for pyscf.ao2mo inside the loaded optimizer module: RUCJ.__init__ calls
    ao2mo.kernel(mf.mol, mo_coeff), which needs AO integrals; a fake-molecule RHF only has mf._eri.
    `kernel` falls back to mf._eri when the molecule has no basis.  Everything else passes through."""

    def __init__(self):
        from pyscf import ao2mo
        self._ao2mo = ao2mo
        self.eri_of_mol = {}

    def __getattr__(self, k):
        return getattr(self._ao2mo, k)

    def kernel(self, mol, mo_coeff, *a, **k):
        eri = self.eri_of_mol.get(id(mol))
        if eri is None:
            return self._ao2mo.kernel(mol, mo_coeff, *a, **k)
        return self._ao2mo.full(eri, mo_coeff)


def load_avijit():
    spec = importlib.util.spec_from_file_location("UCJ_optimizer", AVIJIT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.ao2mo = _Ao2moShim()
    return mod


def make_ucj(mod, mf):
    """mod.UCJ(mf), also for a fake-molecule RHF from `active_mf`."""
    if getattr(mf, "_eri", None) is not None and mf.mol.nbas == 0:
        mod.ao2mo.eri_of_mol[id(mf.mol)] = mf._eri
    return mod.UCJ(mf)


SMALL = {
    "H4_lin_1.5": ("H 0 0 0; H 0 0 1.5; H 0 0 3.0; H 0 0 4.5", "sto-3g"),
    "H4_lin_2.0": ("H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0", "sto-3g"),
    "H4_sq_1.5": ("H 0 0 0; H 1.5 0 0; H 0 1.5 0; H 1.5 1.5 0", "sto-3g"),
    "LiH_1.6": ("Li 0 0 0; H 0 0 1.6", "sto-3g"),
    "LiH_3.0": ("Li 0 0 0; H 0 0 3.0", "sto-3g"),
    "LiH_1.6_631g": ("Li 0 0 0; H 0 0 1.6", "6-31g"),
    "C2_1.25": ("C 0 0 0; C 0 0 1.25", "sto-3g"),
}


def small_mf(name: str, verbose: int = 0):
    from pyscf import gto, scf
    atom, basis = SMALL[name]
    mol = gto.M(atom=atom, basis=basis, unit="Angstrom", symmetry=False, verbose=verbose)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        mf = mf.newton().run()
    return mf


def active_mf(h1: np.ndarray, eri: np.ndarray, nocc: int, const: float, verbose: int = 0):
    """RHF object on an active-space Hamiltonian already expressed in its RHF orbitals.

    Orbitals stay the given ones (mo_coeff = I, mo_occ from aufbau); the SCF is converged in one step
    when the input basis is canonical RHF.  `mol.energy_nuc()` returns `const` (nuclear + frozen core)."""
    from pyscf import ao2mo, gto, scf
    n = h1.shape[0]
    mol = gto.M(verbose=verbose)
    mol.nelectron = 2 * nocc
    mol.nao_nr = lambda *a, **k: n
    mol.incore_anyway = True
    mol.energy_nuc = lambda *a, **k: float(const)
    mf = scf.RHF(mol)
    mf.get_hcore = lambda *a, **k: np.asarray(h1)
    mf.get_ovlp = lambda *a, **k: np.eye(n)
    mf._eri = ao2mo.restore(8, np.asarray(eri), n)
    mf.conv_tol = 1e-11
    dm0 = np.zeros((n, n))
    dm0[np.arange(nocc), np.arange(nocc)] = 2.0
    mf.kernel(dm0)
    return mf


def mo_integrals(mf):
    from pyscf import ao2mo
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    if getattr(mf, "_eri", None) is not None:
        eri = ao2mo.restore(1, ao2mo.full(mf._eri, C), C.shape[1])
    else:
        eri = ao2mo.restore(1, ao2mo.kernel(mf.mol, C), C.shape[1])
    return h1, eri


def ffsim_energy(h1, eri, const, nocc, g, rotation, sign=+1):
    """Exact <psi|H|psi> of exp(i sign X)|HF>, X = rotation-transformed 1/2 sum g n n (Avijit's ansatz)."""
    import ffsim
    n = h1.shape[0]
    ham = ffsim.MolecularHamiltonian(h1, eri, const)
    op = ffsim.UCJOpSpinBalanced(
        diag_coulomb_mats=sign * np.stack([g, g])[None].astype(float),
        orbital_rotations=np.asarray(rotation, dtype=complex).T[None])
    nelec = (nocc, nocc)
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(n, nelec), op, norb=n, nelec=nelec)
    return float(np.vdot(psi, ffsim.linear_operator(ham, norb=n, nelec=nelec) @ psi).real)


def ref_energies(mf, fci: bool = True):
    from pyscf import cc, ci, fci as pfci
    out = dict(e_hf=float(mf.e_tot))
    mycc = cc.CCSD(mf)
    mycc.conv_tol = 1e-9
    mycc.max_cycle = 200
    mycc.kernel()
    out["e_ccsd"] = float(mycc.e_tot)
    out["ccsd_converged"] = bool(mycc.converged)
    myci = ci.CISD(mf)
    myci.kernel()
    out["e_cisd"] = float(myci.e_tot)
    if fci:
        h1, eri = mo_integrals(mf)
        n = h1.shape[0]
        nelec = mf.mol.nelectron
        e, _ = pfci.direct_spin1.kernel(h1, eri, n, nelec, ecore=mf.mol.energy_nuc(), conv_tol=1e-10)
        out["e_fci"] = float(e)
    return out
