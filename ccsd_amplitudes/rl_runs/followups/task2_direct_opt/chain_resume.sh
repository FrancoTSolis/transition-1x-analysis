#!/bin/bash
# Resume (5 Oct, successor agent). After the NOMAD-QSCI run on C2H3N_rxn2858_P has stopped (stop file after 346
# evaluations, then killed; JSON from --finalize): Task-1-protocol scores of the objective-B (QSCI-optimized) parameter sets, then the fixed-size
# subspace scan (task2_dimscan) of both objective-B molecules.  3 SCI threads; scai2 GPU 3.
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
D=rl_runs/followups/task2_direct_opt
until grep -q "queue done" $D/queue_B.log; do sleep 30; done
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=pretrain/.train_venv/bin/python3
echo "=== $(date '+%F %T') scoring objective-B sets"
$PY -m pretrain.followups.task2_score --only qsci --threads 3 --max-mem-gb 2.0 --no-starts \
    --names C2H3N_rxn2858_P C2H3N_rxn2858_TS < /dev/null 2>&1 | grep -v -i warn
echo "=== $(date '+%F %T') dimscan"
$PY -m pretrain.followups.task2_dimscan --names C2H3N_rxn2858_TS C2H3N_rxn2858_P --threads 3 --max-mem-gb 2.0 \
    < /dev/null 2>&1 | grep -v -i warn
echo "=== $(date '+%F %T') chain done"
