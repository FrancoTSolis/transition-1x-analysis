#!/bin/bash
# Task 3: optimize=True labels to convergence (L-BFGS maxiter 5000, square_reg0.005, canonical init) on the 5
# evaluation sets; 8 single-threaded shards on scai2 (CPU budget 8 cores). Small sets first, then n29.
#   CORES="2 7 10 19 23 28 38 41" bash pretrain/followups/task3_run_fits.sh
# Each shard is pinned to one core with taskset: JAX/XLA ignores the thread env vars and otherwise uses 2-4 cores
# per process (the original label runs, run_gen_compressed_scai1.sh, were pinned the same way). Pinning does not
# change the numbers: the 500-iteration snapshots are bit-identical to the existing labels either way.
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
cd $ROOT
CORES=(${CORES:-2 7 10 19 23 28 38 41})
N=${#CORES[@]}
LOG=rl_runs/followups/task3_converged_labels
mkdir -p $LOG
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for i in $(seq 0 $((N-1))); do
  setsid nohup taskset -c ${CORES[$i]} .lucj_venv/bin/python3 -u pretrain/followups/task3_fit_labels.py \
      --names-file pretrain/rl/small_val.txt pretrain/rl/large_val.txt gauge_study/names_norb19.txt \
                   pretrain/rl/n29_val20.txt pretrain/rl/n29_test24.txt \
      --out-dir rhf_targets_compressed_conv --config square_reg0.005 --maxiter 5000 \
      --snap-iters 500 1000 2000 3000 4000 \
      --ref-dirs rhf_targets_compressed_small rhf_targets_compressed_n17 rhf_targets_compressed_n18 \
                 rhf_targets_compressed_n19 rhf_targets_compressed rhf_targets_compressed_n29test \
      --shard $i --n-shards $N >> $LOG/fit_shard$i.log 2>&1 < /dev/null &
done
