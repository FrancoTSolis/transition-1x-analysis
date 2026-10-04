#!/bin/bash
# Debug-partition smoke test of --reward tn in pretrain.rl.grpo_slot (n29, 1 step, tiny groups).
#SBATCH --job-name=mllucj_n29smoke
#SBATCH --account=cla361
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/n29smoke_%j.out
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh; module purge; module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
export OMP_NUM_THREADS=2 TMPDIR=/scratch/$USER/job_$SLURM_JOB_ID; mkdir -p $TMPDIR rl_runs/tn_basis_cache
python -u -m pretrain.rl.grpo_slot --init runs_ot/slotall_T4/best.pt --prefix-T 3 --init-policy rl_runs/grpo_slot_T4p/policy_best.pt \
  --train-names pretrain/rl/n29_smoke_train.txt --val-names pretrain/rl/n29_smoke_val.txt --reward tn --tn-chi 32 \
  --tn-basis-cache rl_runs/tn_basis_cache --group-size 4 --batch-mols 2 --n-workers 12 --worker-threads 1 --device cpu \
  --antithetic --sigma-k 0.01 --sigma-z 0.005 --ppo-epochs 1 --steps 1 --eval-every 1 --out rl_runs/n29_smoke_$SLURM_JOB_ID
echo "done $(date)"
