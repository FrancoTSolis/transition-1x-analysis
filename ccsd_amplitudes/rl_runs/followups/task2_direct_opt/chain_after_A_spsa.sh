#!/bin/bash
# After the A-spsa queue: score the objective-A parameter sets (Task-1 protocol, 1 SCI thread), then start the
# objective-A extension queue.  CPU: 1 core in this slot (B runs hold the other 5).
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
D=rl_runs/followups/task2_direct_opt
until grep -q "queue done" $D/queue_A_spsa.log; do sleep 30; done
echo "=== $(date '+%F %T') scoring objective-A sets"
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
pretrain/.train_venv/bin/python3 -m pretrain.followups.task2_score --only lucj --threads 1 --max-mem-gb 2.0 \
    --names C2H3N_rxn2858_P C2H3N_rxn2858_TS C2H4O_rxn2507_TS < /dev/null 2>&1 | grep -v -i warn
echo "=== $(date '+%F %T') scoring done; starting the extension queue"
bash pretrain/followups/task2_queue2.sh $D/jobs_A_ext.txt > $D/queue_A_ext.log 2>&1 < /dev/null
