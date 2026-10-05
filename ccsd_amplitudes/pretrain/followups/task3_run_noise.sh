#!/bin/bash
# Task 3 extra: 4 more label trajectories at norb 15-17 from the canonical init perturbed by relative 1e-15 Gaussian
# noise (round-off level; seeds 1-4), single-threaded, default L-BFGS-B tolerances, maxiter 5000. Together with the
# stored (multi-thread) and the 1-thread trajectories this gives 6 equally valid labels per molecule.
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
cd $ROOT
LOG=rl_runs/followups/task3_converged_labels
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for s in 1 2 3 4; do
  taskset -c 2 .lucj_venv/bin/python3 -u pretrain/followups/task3_fit_labels.py \
      --names-file pretrain/followups/task3_lists/traj_n15_17.txt --out-dir rhf_targets_compressed_conv/traj_s$s \
      --config square_reg0.005 --maxiter 5000 --snap-iters 500 1000 2000 --x0-noise 1e-15 --x0-seed $s \
      > $LOG/fit_traj_s$s.log 2>&1 < /dev/null
done
