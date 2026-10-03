import os, sys, time, numpy as np
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
os.chdir("/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
from pretrain.rl.hamiltonian import load_hamiltonian
from pretrain.rl.energy import exact_energy, make_ucj_op
name = sys.argv[1]
ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian("rhf_hamiltonians", name)
ini = np.load(f"{os.environ.get('LABELS_DIR', 'rhf_targets_compressed_small')}/square_reg0.005/{name}.npz")
U = ini["U_re"] + 1j * ini["U_im"]; Z = ini["Z"]; t1 = np.load(f"rhf_dataset/{name}.npz")["t1"].astype(float)
t0 = time.time(); E = exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1=t1))
print(f"{name} norb {norb} nelec {nelec}: {time.time()-t0:.1f}s corr% {(e_hf-E)/(e_hf-e_ccsd)*100:.1f}", flush=True)
