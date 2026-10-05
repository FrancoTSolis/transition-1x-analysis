#!/bin/bash
# Task 1 (noiseless QSCI, Lin et al. protocol): PySCF selected-CI diagonalizations of the sampled subspaces on an
# Expanse shared-partition allocation.  Several `task1_sqd_lin diag` workers share the work through O_EXCL claim
# files; each (sample file, protocol) unit is saved after every batch.  The local scai2 workers run disjoint units
# (selected with --protocols / --only-cands / --skip-cands / --seeds), so results never collide when synced back.
#   sync:   python3 expanse/expanse.py put <local samples dir>/ <remote samples dir>/ ; put the two module files
#   submit: sbatch --export=ALL,SPEC="..." pretrain/followups/task1_sqd_diag.slurm
#   SPEC = ';'-separated worker argument groups, e.g.
#          "--protocols lin --only-cands cdf_lin;--protocols n2631g --seeds 0"
#   W workers per group (default 4) x T threads (default 4); cpus-per-task must be >= groups * W * T.
#SBATCH --job-name=t1_sqd_diag
#SBATCH --account=cla361
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/logs/t1_sqd_diag_%j.out
#SBATCH --export=ALL
set -eo pipefail
CODE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
source /etc/profile.d/modules.sh; module purge; module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /home/fsun3/envs/ml_lucj
cd "$CODE"
W=${W:-4}; T=${T:-4}
export OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T RAYON_NUM_THREADS=$T NUMBA_NUM_THREADS=$T
export JAX_PLATFORMS=cpu XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"
LOG=rl_runs/followups/task1_sqd_vs_lin/expanse_$SLURM_JOB_ID
mkdir -p $LOG
echo "node $(hostname) job $SLURM_JOB_ID cpus $SLURM_CPUS_PER_TASK  W=$W T=$T  SPEC=$SPEC  $(date)"
IFS=';' read -ra GROUPS_ <<< "${SPEC:---protocols lin}"
pids=()
g=0
for args in "${GROUPS_[@]}"; do
  for i in $(seq 1 $W); do
    python -u -m pretrain.followups.task1_sqd_lin diag --threads $T $args --idle-exit-s ${IDLE:-900} \
      --stop-file $LOG/STOP > $LOG/g${g}_w$i.log 2>&1 &
    pids+=($!)
  done
  g=$((g + 1))
done
wait "${pids[@]}"
echo "done $(date)"
