#!/bin/bash
# GRPO on n29 molecules (norb 29, 58 qubits) with the tensor-network (MPS) LUCJ energy as reward
# (pretrain/rl/tn_energy.py: split-localized basis, exact square-mask factors, block2 energy), one Expanse node:
# policy driver on CPU + 96 single-threaded TN reward workers (~2.7 min per energy at chi 64, ~1.6-2 GB each with
# one cached engine; more workers or cached engines ran out of the node's 240 GB).
#   sbatch --export=ALL,RUN=grpo_n29_tn,INITPOL=rl_runs/grpo_large_rl4/policy_snap_for_n29.pt expanse/grpo_n29_tn.slurm
# (TNIMPL=v1, the default here, is the frozen zip-up engine of the Oct-3 run; TNIMPL=current uses tn_energy.py.)
#SBATCH --job-name=mllucj_n29tn
#SBATCH --account=cla361
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=06:00:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/n29tn_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh; module purge; module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
# everything single-threaded: spawned reward workers inherit this before they can pin themselves (numpy/MKL read
# it at import); the driver gets its torch threads through --driver-threads
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export TMPDIR=/scratch/$USER/job_$SLURM_JOB_ID; mkdir -p $TMPDIR      # node-local block2 scratch
RUN=${RUN:-grpo_n29_tn_$SLURM_JOB_ID}
mkdir -p rl_runs/tn_basis_cache
echo "node $(hostname) job $SLURM_JOB_ID run $RUN $(date)"
python -u -m pretrain.rl.grpo_slot --init runs_ot/slotall_T4/best.pt --prefix-T 3 --init-policy ${INITPOL:-rl_runs/grpo_slot_T4p/policy_best.pt} \
  --train-names pretrain/rl/n29_train40.txt --val-names pretrain/rl/n29_val20.txt \
  --reward tn --tn-chi ${CHI:-64} --tn-impl ${TNIMPL:-v1} --tn-basis-cache rl_runs/tn_basis_cache --tn-stack-mem-gb 1 \
  --group-size 8 --batch-mols 12 --n-workers 96 --worker-threads 1 --device cpu --driver-threads 8 --tn-cache-items 1 \
  --sigma-k 0.01 --sigma-z 0.005 --antithetic --ppo-epochs 1 --lr ${LR:-1e-5} --steps 30 --eval-every 5 \
  --out rl_runs/$RUN ${ARGS:-}
echo "done $(date)"
