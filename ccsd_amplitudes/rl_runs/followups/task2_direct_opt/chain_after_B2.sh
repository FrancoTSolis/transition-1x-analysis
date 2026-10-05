#!/bin/bash
# After the second QSCI NOMAD run: SPSA on the QSCI objective (2 SCI threads, 3 h cap each).
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
D=rl_runs/followups/task2_direct_opt
until grep -q "queue done" $D/queue_B2.log; do sleep 30; done
bash pretrain/followups/task2_queue2.sh $D/jobs_B3.txt > $D/queue_B3.log 2>&1 < /dev/null
