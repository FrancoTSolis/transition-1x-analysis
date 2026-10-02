import os, sys, collections
os.environ["OMP_NUM_THREADS"] = "2"
import numpy as np
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
import ffsim, quimb.tensor as qtn
from qiskit.circuit import QuantumCircuit, QuantumRegister
from qiskit_quimb import quimb_circuit
from qiskit.primitives import BitArray
from qiskit_addon_sqd.fermion import solve_fermion, diagonalize_fermionic_hamiltonian
from pretrain.rl.hamiltonian import load_hamiltonian
from pretrain.rl.energy import interaction_pairs, lucj_circuit
name = "C2H3N_rxn2857_TS"
ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian("rhf_hamiltonians", name)
d = np.load(f"rhf_dataset/{name}.npz"); t2 = d["t2"].astype(float); t1 = d["t1"].astype(float)
# 1. ffsim's own HF bitstring convention
hf = ffsim.hartree_fock_state(norb, nelec)
s_ffsim = ffsim.sample_state_vector(hf, norb=norb, nelec=nelec, shots=3, seed=0, bitstring_type=ffsim.BitstringType.STRING)
print("ffsim HF string:", s_ffsim[0], "len", len(s_ffsim[0]))
# 2. quimb sample of the HF-only circuit
q = QuantumRegister(2*norb); c = QuantumCircuit(q); c.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), q)
qc = quimb_circuit(c.decompose(reps=2), quimb_circuit_class=qtn.CircuitMPS, max_bond=8)
raw = list(qc.sample(5, seed=0))
print("quimb raw HF sample:", raw[0]); print("reversed          :", raw[0][::-1])
# 3. exact-DF LUCJ circuit, chi=64
op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, t1=t1, n_reps=2, interaction_pairs=interaction_pairs("square", norb))
qc = quimb_circuit(lucj_circuit(op, norb, nelec), quimb_circuit_class=qtn.CircuitMPS, max_bond=64)
raw = list(qc.sample(2000, seed=0))
for label, ss in [("reversed", [s[::-1] for s in raw]), ("raw", raw)]:
    cnt = collections.Counter(ss)
    top = cnt.most_common(3)
    pn = [(s[:norb].count("1"), s[norb:].count("1")) for s in ss]
    ok = sum(1 for a, b in pn if (a, b) == nelec) / len(ss)
    print(f"[{label}] top: {top}  frac with (n_beta,n_alpha)=={nelec}: {ok:.3f}")
    ba = BitArray.from_samples(ss, num_bits=2*norb)
    try:
        res = diagonalize_fermionic_hamiltonian(ham.one_body_tensor.real, ham.two_body_tensor.real, ba, samples_per_batch=300,
                norb=norb, nelec=nelec, num_batches=1, max_iterations=1, symmetrize_spin=True, seed=0)
        print(f"   SQD it=1: E={res.energy+ham.constant:.5f} corr%={(e_hf-(res.energy+ham.constant))/(e_hf-e_ccsd)*100:.1f} dim={res.sci_state.amplitudes.shape}")
    except Exception as ex:
        print("   SQD failed:", ex)
# 4. exact statevector samples from ffsim for the same op (ground truth distribution)
psi = ffsim.apply_unitary(hf, op, norb=norb, nelec=nelec)
ss = ffsim.sample_state_vector(psi, norb=norb, nelec=nelec, shots=2000, seed=0, bitstring_type=ffsim.BitstringType.STRING)
print("ffsim exact top:", collections.Counter(ss).most_common(3))
ba = BitArray.from_samples(ss, num_bits=2*norb)
res = diagonalize_fermionic_hamiltonian(ham.one_body_tensor.real, ham.two_body_tensor.real, ba, samples_per_batch=300,
        norb=norb, nelec=nelec, num_batches=1, max_iterations=1, symmetrize_spin=True, seed=0)
print(f"   SQD on exact samples it=1: E={res.energy+ham.constant:.5f} corr%={(e_hf-(res.energy+ham.constant))/(e_hf-e_ccsd)*100:.1f} dim={res.sci_state.amplitudes.shape}")
