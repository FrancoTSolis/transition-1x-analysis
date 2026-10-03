#!/bin/bash
# Variant: policy = 3 frozen recycles of the 4-recycle model + trainable final recycle (pretrain.rl.grpo_slot --prefix-T 3).
# GRPO energy fine-tuning of the chemistry-frame slot model on Expanse (one node; policy on CPU, reward pool of
# single-threaded workers).  Exact statevector rewards are memory-bound (norb 16: ~10-12 GB per energy), so the
# worker count follows memory, not cores; MPS+SQD rewards (--reward mps_sqd) need far less memory.
#   sbatch --export=ALL,RUN=grpo_slot_v1,ARGS="--steps 40" expanse/grpo_slot.slurm
#SBATCH --job-name=mllucj_grpo
#SBATCH --account=cla361
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=05:00:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/grpo_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh
module purge
module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4   # policy (torch CPU); workers pin themselves
RUN=${RUN:-grpo_slot_$SLURM_JOB_ID}
echo "node $(hostname) job $SLURM_JOB_ID run $RUN $(date)"
python -u -m pretrain.rl.grpo_slot --init runs_ot/slotall_T4/best.pt --prefix-T 3 \
  --train-names pretrain/rl/small_train.txt --val-names pretrain/rl/small_val.txt \
  --group-size 8 --batch-mols 4 --n-workers 16 --worker-threads 7 --device cpu \
  --sigma-k 0.01 --sigma-z 0.005 --antithetic --steps 30 --eval-every 5 \
  --out rl_runs/$RUN ${ARGS:-}
echo "done $(date)"
