#!/bin/bash
# launch_slot.sh <gpu> <run_name> [train_slot args]   (local scai2 only)
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
g=$1; name=$2; shift 2
CUDA_VISIBLE_DEVICES=$g OMP_NUM_THREADS=3 setsid nohup pretrain/.train_venv/bin/python3 -u -m pretrain.opt_true.train_slot \
  --out runs_ot/$name "$@" > runs_ot/$name.log 2>&1 < /dev/null &
echo "launched $name on GPU $g (pid $!)"
