#!/bin/bash
# 2-step GRPO check on the debug partition with the instant residual reward (verifies the Expanse setup).
#SBATCH --job-name=mllucj_grpotest
#SBATCH --account=cla361
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/grpotest_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh
module purge
module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
export OMP_NUM_THREADS=4
python -u -m pretrain.rl.grpo_slot --init runs_ot/slotall_T1/best.pt --train-names pretrain/rl/small_train.txt \
  --val-names pretrain/rl/small_val.txt --reward resid --group-size 8 --batch-mols 2 --steps 2 --eval-every 2 \
  --n-workers 1 --device cpu --antithetic --sigma-k 0.01 --sigma-z 0.005 --out rl_runs/grpotest_$SLURM_JOB_ID
echo "done $(date)"
