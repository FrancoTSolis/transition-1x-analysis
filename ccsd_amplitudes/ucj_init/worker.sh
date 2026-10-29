#!/bin/bash
# usage: worker.sh <jobs file> <gpu index>   (run from any scai host; the tree is on shared storage)
# Each line of the jobs file is the argument list of ucj_init.opt_exact.  A job is claimed by an atomic mkdir.
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
JOBS=$1; GPU=$2
LOCK=rl_runs/ucj_init/exact_opt/locks; mkdir -p $LOCK
export CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
while read -r line; do
  [ -z "$line" ] && continue
  key=$(echo "$line" | md5sum | cut -c1-12)
  mkdir $LOCK/$key 2>/dev/null || continue
  echo "$(date +%T) $(hostname) gpu$GPU start: $line"
  pretrain/.train_venv/bin/python -m ucj_init.opt_exact $line 2>&1 | grep -v -i warn | tail -1
  echo "$(date +%T) done: $line"
done < $JOBS
