#!/bin/bash
# Debug-partition cost test of the TN reward at n29: 8 single-threaded energies in parallel (chi 64).
#SBATCH --job-name=mllucj_tnbench
#SBATCH --account=cla361
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/tnbench_%j.out
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh; module purge; module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 TMPDIR=/scratch/$USER/job_$SLURM_JOB_ID
mkdir -p $TMPDIR
echo "node $(hostname) $(date)"
for n in $(head -8 pretrain/rl/n29_val20.txt); do python -u expanse/tn_bench.py --name $n --chi ${CHI:-64} --out rl_runs/tn_bench_$SLURM_JOB_ID.jsonl & done
wait
echo "done $(date)"
