#!/bin/bash
# Exact energies for larger active spaces (norb 17: ~42 GB and ~2 core-h per energy): few multi-threaded workers.
#   sbatch --export=ALL,DUMP=runs_ot/energy_tasks/n17.pkl,OUT=pretrain/opt_true/results/energy_n17.json expanse/energy_dump_big.slurm
#SBATCH --job-name=mllucj_energy17
#SBATCH --account=cla361
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/energy17_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh
module purge
module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
export OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= ENERGY_THREADS=${ENERGY_THREADS:-25}
echo "node $(hostname) job $SLURM_JOB_ID $(date)"
python -u -m pretrain.opt_true.eval_slot_energy --from-dump "$DUMP" --n-procs ${NPROCS:-5} --out "$OUT" --skip-keys ${SKIP:-}
echo "done $(date)"
