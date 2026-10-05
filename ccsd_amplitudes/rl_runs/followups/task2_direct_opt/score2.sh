#!/bin/bash
# Resume (5 Oct, successor agent): Task-1-protocol scores of the objective-A sets of the six extension molecules and
# the long label run (2 SCI threads; scai2 GPU 3).
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
echo "=== $(date '+%F %T') scoring ($*)"
pretrain/.train_venv/bin/python3 -m pretrain.followups.task2_score "$@" < /dev/null 2>&1 | grep -v -i warn
echo "=== $(date '+%F %T') scoring done"
