import os, time, sys
os.environ["OMP_NUM_THREADS"] = "4"
import numpy as np
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
import ffsim
from pretrain.rl.hamiltonian import load_hamiltonian
from pretrain.rl.energy import mps_sample, sqd_energy_from_samples, interaction_pairs, exact_energy
name = sys.argv[1]; conn = "square"
ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian("rhf_hamiltonians", name)
d = np.load(f"rhf_dataset/{name}.npz"); t2 = d["t2"].astype(float); t1 = d["t1"].astype(float)
op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, t1=t1, n_reps=2, interaction_pairs=interaction_pairs(conn, norb))
print(f"{name} norb={norb} nelec={nelec} corr={e_ccsd-e_hf:.4f}", flush=True)
for chi, shots in [(16, 500), (32, 1000), (64, 1000), (64, 3000)]:
    t0 = time.time(); samples, info = mps_sample(op, norb, nelec, shots=shots, max_bond=chi, seed=0)
    print(f"  MPS chi={chi} shots={shots}: t_mps={info['t_mps']:.1f}s t_sample={info['t_sample']:.1f}s uniq={len(set(samples))}", flush=True)
    for (nb, it, spb) in [(1, 1, 200), (1, 2, 300), (3, 5, 300)]:
        t0 = time.time()
        try:
            e, res = sqd_energy_from_samples(ham, norb, nelec, samples, samples_per_batch=spb, n_batches=nb, max_iterations=it, seed=0)
            print(f"    SQD nb={nb} it={it} spb={spb}: E={e:.6f} corr%={(e_hf-e)/(e_hf-e_ccsd)*100:6.1f} dim={res.sci_state.amplitudes.shape} ({time.time()-t0:.1f}s)", flush=True)
        except Exception as ex:
            print("    SQD failed", type(ex).__name__, ex, flush=True)
