import os, sys, time, collections
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
import ffsim, quimb.tensor as qtn
from qiskit_quimb import quimb_circuit
from pretrain.rl.hamiltonian import load_hamiltonian
from pretrain.rl.energy import interaction_pairs, lucj_circuit
name = sys.argv[1] if len(sys.argv) > 1 else "C2H3N_rxn2857_TS"
ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian("rhf_hamiltonians", name)
d = np.load(f"rhf_dataset/{name}.npz"); t2 = d["t2"].astype(float); t1 = d["t1"].astype(float)
op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, t1=t1, n_reps=2, interaction_pairs=interaction_pairs("square", norb))
circ = lucj_circuit(op, norb, nelec)
print(name, "norb", norb, "qubits", 2*norb, "gates:", collections.Counter(g.operation.name for g in circ.data).most_common(8), flush=True)
hf_str = ffsim.sample_state_vector(ffsim.hartree_fock_state(norb, nelec), norb=norb, nelec=nelec, shots=1, seed=0, bitstring_type=ffsim.BitstringType.STRING)[0]
if norb <= 16:
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, nelec), op, norb=norb, nelec=nelec)
    p_hf_exact = abs(psi[ffsim.dim(norb, nelec)-1 if False else 0])**2  # placeholder
    ss = ffsim.sample_state_vector(psi, norb=norb, nelec=nelec, shots=4000, seed=0, bitstring_type=ffsim.BitstringType.STRING)
    print(f"exact: P(HF) ~ {collections.Counter(ss)[hf_str]/4000:.3f}", flush=True)
for chi, cutoff in [(64, 1e-8), (128, 1e-8), (256, 1e-10), (None, 1e-12)]:
    t0 = time.time()
    try:
        qc = quimb_circuit(circ, quimb_circuit_class=qtn.CircuitMPS, max_bond=chi, cutoff=cutoff, progbar=False)
        t1_ = time.time()
        psi_mps = qc.psi
        nrm = float(abs(psi_mps.norm()))
        amp = qc.amplitude(hf_str[::-1])   # quimb order = reversed qiskit string
        t2_ = time.time()
        smp = [s[::-1] for s in qc.sample(1000, seed=0)]
        frac_hf = sum(1 for s in smp if s == hf_str) / len(smp)
        frac_pn = sum(1 for s in smp if (s[:norb].count("1"), s[norb:].count("1")) == nelec) / len(smp)
        print(f"chi={chi} cutoff={cutoff}: maxbond={psi_mps.max_bond()} norm={nrm:.4f} |<HF|psi>|^2={abs(amp)**2:.4f} "
              f"sampled P(HF)={frac_hf:.3f} P(particle-number ok)={frac_pn:.3f}  t_build={t1_-t0:.1f}s t_sample={time.time()-t2_:.1f}s", flush=True)
    except Exception as ex:
        print(f"chi={chi}: FAILED {type(ex).__name__}: {ex}", flush=True)
