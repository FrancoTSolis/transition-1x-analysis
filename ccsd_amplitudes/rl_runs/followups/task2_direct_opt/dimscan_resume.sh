#!/bin/bash
# Resume (5 Oct, successor agent): fixed-size subspace scan (task2_dimscan, k = 150..800 most probable strings) of the
# two objective-B molecules, after the objective-A scoring (score2.sh) has freed its 2 cores.  scai2 GPU 3.
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
D=rl_runs/followups/task2_direct_opt
until grep -q "scoring done" $D/queue_score2.log; do sleep 30; done
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
echo "=== $(date '+%F %T') dimscan"
pretrain/.train_venv/bin/python3 -u -m pretrain.followups.task2_dimscan --names C2H3N_rxn2858_TS C2H3N_rxn2858_P \
    --threads 2 --max-mem-gb 2.0 < /dev/null 2>&1 | grep --line-buffered -v -i warn
echo "=== $(date '+%F %T') dimscan done"
