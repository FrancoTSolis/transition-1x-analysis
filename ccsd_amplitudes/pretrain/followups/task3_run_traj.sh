#!/bin/bash
# Task 3 extra: label run-to-run spread. The same fits (canonical init, square_reg0.005, default L-BFGS-B tolerances,
# maxiter 5000) under different XLA CPU thread-pool sizes (= CPU affinity size), which changes only the float
# summation order. The stored labels at norb 15-17 match the unpinned (48-core pool) run bit-for-bit.
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
cd $ROOT
LOG=rl_runs/followups/task3_converged_labels
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for spec in "p1:2" "p2:2,7"; do
  tag=${spec%%:*}; cores=${spec#*:}
  taskset -c $cores .lucj_venv/bin/python3 -u pretrain/followups/task3_fit_labels.py \
      --names-file pretrain/followups/task3_lists/traj_n15_17.txt --out-dir rhf_targets_compressed_conv/traj_$tag \
      --config square_reg0.005 --maxiter 5000 --snap-iters 500 1000 2000 \
      --ref-dirs rhf_targets_compressed_small rhf_targets_compressed_n17 > $LOG/fit_traj_$tag.log 2>&1 < /dev/null
done
