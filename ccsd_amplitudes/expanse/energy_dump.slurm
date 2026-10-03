#!/bin/bash
# Exact LUCJ energies for pickled candidate tasks (pretrain.opt_true.eval_slot_energy --dump) on Expanse.
# Memory-bound (norb 16 ~10-12 GB per energy): 20 single-threaded workers in 240 GB.
#   sbatch --export=ALL,DUMP=runs_ot/energy_tasks/all1.pkl,OUT=pretrain/opt_true/results/energy_small_all1.json,SKIP=label expanse/energy_dump.slurm
#SBATCH --job-name=mllucj_energy
#SBATCH --account=cla361
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=240G
#SBATCH --time=02:00:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/energy_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh
module purge
module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=
echo "node $(hostname) job $SLURM_JOB_ID $(date)"
python -u -m pretrain.opt_true.eval_slot_energy --from-dump "$DUMP" --n-procs ${NPROCS:-20} --out "$OUT" --skip-keys ${SKIP:-}
echo "done $(date)"
