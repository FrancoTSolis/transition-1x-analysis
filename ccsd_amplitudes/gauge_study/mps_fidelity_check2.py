"""MPS fidelity of the LUCJ circuit vs qubit ordering / MPS class.
Default JW order puts alpha orbitals on qubits 0..n-1 and beta on n..2n-1, so the
alpha-beta diagonal-Coulomb gates span n qubits (long-range swaps + truncation).
Test: (a) interleaved order (alpha_p, beta_p adjacent), (b) CircuitPermMPS."""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
import ffsim, quimb.tensor as qtn
from qiskit.circuit import QuantumCircuit
from qiskit_quimb import quimb_circuit
from pretrain.rl.hamiltonian import load_hamiltonian
from pretrain.rl.energy import interaction_pairs, lucj_circuit
name = sys.argv[1] if len(sys.argv) > 1 else "C2H3N_rxn2857_TS"
ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian("rhf_hamiltonians", name)
d = np.load(f"rhf_dataset/{name}.npz"); t2 = d["t2"].astype(float); t1 = d["t1"].astype(float)
op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, t1=t1, n_reps=2, interaction_pairs=interaction_pairs("square", norb))
circ = lucj_circuit(op, norb, nelec)
hf_q = ffsim.sample_state_vector(ffsim.hartree_fock_state(norb, nelec), norb=norb, nelec=nelec, shots=1, seed=0,
                                 bitstring_type=ffsim.BitstringType.STRING)[0]   # qiskit order (qubit 2n-1 ... 0)
hf_bits = hf_q[::-1]                                                              # index by qubit
# interleaved layout: logical qubit q -> physical position perm[q]; alpha_p -> 2p, beta_p -> 2p+1
perm = [2 * q if q < norb else 2 * (q - norb) + 1 for q in range(2 * norb)]
circ_il = QuantumCircuit(2 * norb); circ_il.compose(circ, qubits=perm, inplace=True)
hf_il = [None] * (2 * norb)
for q in range(2 * norb):
    hf_il[perm[q]] = hf_bits[q]
hf_il = "".join(hf_il)
print(name, "norb", norb, flush=True)
for label, c, hfs, cls, chis in [
        ("default  CircuitMPS", circ, hf_bits, qtn.CircuitMPS, [64, 128]),
        ("interleaved CircuitMPS", circ_il, hf_il, qtn.CircuitMPS, [64, 128, 256]),
        ("default  CircuitPermMPS", circ, hf_bits, qtn.CircuitPermMPS, [64, 128]),
]:
    for chi in chis:
        t0 = time.time()
        try:
            qc = quimb_circuit(c, quimb_circuit_class=cls, max_bond=chi, cutoff=1e-10, progbar=False)
            tb = time.time() - t0
            nrm = float(abs(qc.psi.norm()))
            amp = qc.amplitude(hfs) if cls is qtn.CircuitMPS else qc.amplitude(hfs)
            print(f"  {label:24s} chi={chi:4d}: norm={nrm:.4f}  |<HF|psi>|^2={abs(amp)**2:.4f} (exact 0.968)  "
                  f"maxbond={qc.psi.max_bond()}  t_build={tb:.0f}s", flush=True)
        except Exception as ex:
            print(f"  {label:24s} chi={chi}: FAILED {type(ex).__name__}: {ex}", flush=True)
