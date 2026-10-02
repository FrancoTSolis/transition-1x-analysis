"""Energy evaluation of LUCJ (U, Z) parameters -- the RL reward.

Two evaluators, both returning the total electronic energy in Hartree of the
state  prod_k U_k exp(i J_k) U_k^dag  [x final orbital rotation from t1]  |HF>:

  exact_energy(...)      : ffsim statevector, feasible for norb <= ~16
                           (dim = C(norb, nocc)^2; norb=15/nocc=8 -> 4e7).
  mps_sqd_energy(...)    : Jordan-Wigner circuit -> quimb MPS (max_bond,
                           cutoff) -> sample -> SQD (qiskit-addon-sqd,
                           pyscf selected-CI solver).  This is the reward of
                           kevinsung/lucj's lucj_sqd_quimb_task_nomad.py,
                           minus the NOMAD black-box optimizer around it.

(U, Z) -> ffsim.UCJOpSpinBalanced is built with the same convention as the
label pipeline: diag_coulomb_mats = stack([Z, Z], axis=1) masked to the
alpha-alpha / alpha-beta interaction pairs of the chosen connectivity.
"""
from __future__ import annotations

import time

import numpy as np


def interaction_pairs(connectivity: str, norb: int):
    if connectivity == "all-to-all":
        return None, None
    if connectivity == "square":
        return [(p, p + 1) for p in range(norb - 1)], [(p, p) for p in range(norb)]
    if connectivity == "hex":
        return ([(p, p + 1) for p in range(norb - 1)],
                [(p, p) for p in range(norb) if p % 2 == 0])
    if connectivity == "heavy-hex":
        return ([(p, p + 1) for p in range(norb - 1)],
                [(p, p) for p in range(norb) if p % 4 == 0])
    raise ValueError(connectivity)


def make_ucj_op(Z: np.ndarray, U: np.ndarray, connectivity: str,
                t1: np.ndarray | None = None):
    """(Z, U) stacks -> ffsim.UCJOpSpinBalanced with the LUCJ sparsity applied."""
    import ffsim
    from ffsim.variational.util import orbital_rotation_from_t1_amplitudes

    n_reps, norb, _ = Z.shape
    pairs_aa, pairs_ab = interaction_pairs(connectivity, norb)
    J = np.stack([Z, Z], axis=1).astype(float)
    for spin, pairs in ((0, pairs_aa), (1, pairs_ab)):
        if pairs is not None:
            mask = np.zeros((norb, norb), dtype=bool)
            if pairs:
                r, c = zip(*pairs)
                mask[list(r), list(c)] = True
                mask[list(c), list(r)] = True
            J[:, spin] *= mask
    final = orbital_rotation_from_t1_amplitudes(t1) if t1 is not None else None
    return ffsim.UCJOpSpinBalanced(
        diag_coulomb_mats=J, orbital_rotations=U.astype(complex),
        final_orbital_rotation=final)


def z_from_op(op, connectivity: str) -> np.ndarray:
    """Inverse of make_ucj_op's masking: recombine ffsim's per-spin diagonal
    Coulomb blocks into the single Z used by the DF label pipeline."""
    J = op.diag_coulomb_mats
    if connectivity == "all-to-all":
        return J[:, 0].copy()
    norb = J.shape[-1]
    pairs_aa, pairs_ab = interaction_pairs(connectivity, norb)
    Z = np.zeros_like(J[:, 0])
    for spin, pairs in ((0, pairs_aa), (1, pairs_ab)):
        r, c = zip(*pairs)
        Z[:, list(r), list(c)] = J[:, spin][:, list(r), list(c)]
        Z[:, list(c), list(r)] = J[:, spin][:, list(c), list(r)]
    return Z


# ----------------------------------------------------------------- exact

def exact_energy(ham, norb: int, nelec: tuple[int, int], op) -> float:
    import ffsim
    ref = ffsim.hartree_fock_state(norb, nelec)
    psi = ffsim.apply_unitary(ref, op, norb=norb, nelec=nelec)
    linop = ffsim.linear_operator(ham, norb=norb, nelec=nelec)
    return float(np.vdot(psi, linop @ psi).real)


# --------------------------------------------------------------- MPS + SQD

def lucj_circuit(op, norb: int, nelec: tuple[int, int]):
    import ffsim
    from qiskit.circuit import QuantumCircuit, QuantumRegister
    qubits = QuantumRegister(2 * norb)
    circ = QuantumCircuit(qubits)
    circ.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qubits)
    circ.append(ffsim.qiskit.UCJOpSpinBalancedJW(op), qubits)
    return circ.decompose(reps=2)


def mps_sample(op, norb: int, nelec: tuple[int, int], *, shots: int,
               max_bond: int, cutoff: float = 1e-8, seed: int = 0,
               perm_mps: bool = False):
    """Sample bitstrings (little-endian ffsim convention) from the MPS-simulated
    LUCJ circuit.  Returns (list_of_bitstrings, dict_timing)."""
    import quimb.tensor as qtn
    from qiskit_quimb import quimb_circuit

    t0 = time.time()
    circ = lucj_circuit(op, norb, nelec)
    qc = quimb_circuit(
        circ,
        quimb_circuit_class=qtn.CircuitPermMPS if perm_mps else qtn.CircuitMPS,
        max_bond=max_bond, cutoff=cutoff, progbar=False)
    t1 = time.time()
    samples = [s[::-1] for s in qc.sample(shots, seed=seed)]
    t2 = time.time()
    info = dict(t_mps=t1 - t0, t_sample=t2 - t1,
                max_bond_reached=int(qc.psi.max_bond()))
    return samples, info


def sqd_energy_from_samples(ham, norb: int, nelec: tuple[int, int], samples,
                            *, samples_per_batch: int = 300, n_batches: int = 3,
                            max_iterations: int = 5, max_dim: int | None = None,
                            symmetrize_spin: bool = True, seed: int = 0,
                            energy_tol: float = 1e-6):
    """SQD (self-consistent configuration recovery) energy from bitstrings."""
    from qiskit.primitives import BitArray
    from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian

    bit_array = BitArray.from_samples(samples, num_bits=2 * norb)
    rng = np.random.default_rng(seed)
    result = diagonalize_fermionic_hamiltonian(
        ham.one_body_tensor.real, ham.two_body_tensor.real, bit_array,
        samples_per_batch=samples_per_batch, norb=norb, nelec=nelec,
        num_batches=n_batches, energy_tol=energy_tol,
        max_iterations=max_iterations, symmetrize_spin=symmetrize_spin,
        max_dim=max_dim, seed=rng)
    return float(result.energy + ham.constant), result


def mps_sqd_energy(ham, norb, nelec, op, *, shots=2000, max_bond=64,
                   cutoff=1e-8, seed=0, **sqd_kw):
    samples, info = mps_sample(op, norb, nelec, shots=shots, max_bond=max_bond,
                               cutoff=cutoff, seed=seed)
    t0 = time.time()
    e, _ = sqd_energy_from_samples(ham, norb, nelec, samples, seed=seed, **sqd_kw)
    info["t_sqd"] = time.time() - t0
    info["n_unique"] = len(set(samples))
    return e, info
