#!/bin/bash
# Task 1: PySCF FCI references (direct_spin0, conv_tol 1e-10) for given molecules, one process per molecule.
# Each process writes its own JSON (fci_<name>.json) so parallel runs never race on one file; merge locally.
#   sbatch --export=ALL,NAMES="C3H4_rxn2391_P C3H4_rxn2391_TS",T=16 pretrain/followups/task1_fci.slurm
#SBATCH --job-name=t1_fci
#SBATCH --account=cla361
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/t1_fci_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh; module purge; module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
T=${T:-16}
export OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T RAYON_NUM_THREADS=$T NUMBA_NUM_THREADS=$T
export JAX_PLATFORMS=cpu
OUT=pretrain/opt_true/results/followups/task1_sqd_vs_lin/expanse_fci
LOG=rl_runs/followups/task1_sqd_vs_lin/expanse_fci_$SLURM_JOB_ID
mkdir -p $OUT $LOG
echo "node $(hostname) job $SLURM_JOB_ID NAMES=$NAMES T=$T $(date)"
pids=()
for n in $NAMES; do
  python -u - "$n" "$T" "$OUT" > $LOG/$n.log 2>&1 <<'EOF' &
import json, socket, sys, time
import numpy as np
sys.path.insert(0, ".")
from pretrain.followups.task1_sqd_lin import load_ham
from pyscf import fci, lib
name, T, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
lib.num_threads(T)
h = load_ham(name)
s = fci.direct_spin0.FCI()
s.max_memory, s.conv_tol, s.max_cycle, s.verbose = 30000, 1e-10, 200, 5
t0 = time.time()
e, c = s.kernel(h["h1"], h["h2"], h["norb"], h["nelec"], ecore=h["const"])
ss, _ = s.spin_square(c, h["norb"], h["nelec"])
rec = {name: dict(e_fci=float(e), s2=float(ss), conv_tol=1e-10, t_s=time.time() - t0, threads=T,
                  solver="pyscf.fci.direct_spin0", converged=bool(np.all(s.converged)), dim=int(c.size),
                  c_hf=float(abs(c.reshape(-1)[0])), host=socket.gethostname())}
json.dump(rec, open(f"{out}/fci_{name}.json", "w"), indent=1)
print(rec)
EOF
  pids+=($!)
done
wait "${pids[@]}"
echo "done $(date)"
