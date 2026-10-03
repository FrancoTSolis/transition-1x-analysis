#!/bin/bash
# Debug-partition benchmark of energy evaluators on Expanse (cheap: <= 32 cores x 30 min).
#   sbatch expanse/bench_energy.slurm
#SBATCH --job-name=mllucj_bench
#SBATCH --account=cla361
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=00:30:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/bench_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh
module purge
module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
echo "node $(hostname) cores $SLURM_CPUS_PER_TASK $(date)"
python -u - <<'EOF'
import os, time, numpy as np, sys
sys.path.insert(0, ".")
import ffsim
from pretrain.rl.hamiltonian import load_hamiltonian
from pretrain.rl.energy import exact_energy, make_ucj_op, interaction_pairs, mps_sample, sqd_energy_from_samples
name = [l.strip() for l in open("gauge_study/names_norb15.txt")][2]
ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian("rhf_hamiltonians", name)
d = np.load(f"rhf_dataset/{name}.npz"); t2 = d["t2"].astype(float); t1 = d["t1"].astype(float)
op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, t1=t1, n_reps=2, interaction_pairs=interaction_pairs("square", norb))
for nt in (1, 4, 8, 16, 32):
    os.environ["OMP_NUM_THREADS"] = str(nt)
    try:
        import threadpoolctl
        ctx = threadpoolctl.threadpool_limits(nt)
    except Exception:
        ctx = None
    t0 = time.time(); E = exact_energy(ham, norb, nelec, op)
    print(f"exact  threads={nt:2d}: {time.time()-t0:6.1f}s  corr% {(e_hf-E)/(e_hf-e_ccsd)*100:.1f}", flush=True)
for chi in (128, None):
    t0 = time.time()
    s, info = mps_sample(op, norb, nelec, shots=2000, max_bond=chi or 4096, cutoff=1e-10, seed=0)
    t1_ = time.time()
    e, _ = sqd_energy_from_samples(ham, norb, nelec, s, samples_per_batch=300, n_batches=3, max_iterations=5, seed=0)
    valid = sum(1 for x in s if (x[:norb].count("1"), x[norb:].count("1")) == nelec) / len(s)
    print(f"MPS chi={chi}: build+sample {t1_-t0:6.1f}s  maxbond {info['max_bond_reached']}  valid {valid:.2f}  SQD {time.time()-t1_:5.1f}s  corr% {(e_hf-e)/(e_hf-e_ccsd)*100:.1f}", flush=True)
EOF
echo "done $(date)"
