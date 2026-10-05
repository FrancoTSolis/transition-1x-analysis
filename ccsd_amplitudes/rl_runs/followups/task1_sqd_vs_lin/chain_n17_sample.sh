#!/bin/bash
# wait for the small-val sampler, then sample the 12 norb-17 C2HNO molecules on the same GPU (scai2 GPU 2)
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
while kill -0 388191 2>/dev/null; do sleep 30; done
export CUDA_VISIBLE_DEVICES=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 RAYON_NUM_THREADS=4 NUMBA_NUM_THREADS=4
exec pretrain/.train_venv/bin/python3 -u -m pretrain.followups.task1_sqd_lin sample \
  --cands runs_ot/energy_tasks/task1_sqd_vs_lin/cands_n17c2hno.pkl --seeds 0 1 2 3 4 --threads 4 --ci-workers 5 --max-mem-gb 10
