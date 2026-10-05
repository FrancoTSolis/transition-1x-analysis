#!/bin/bash
# push new sample files to Expanse and pull results every 10 min, until STOP_sync exists (max 8 h)
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
L=rl_runs/followups/task1_sqd_vs_lin
end=$(( $(date +%s) + 8*3600 ))
while [ ! -e $L/STOP_sync ] && [ $(date +%s) -lt $end ]; do
  echo "== $(date)"; timeout 900 bash pretrain/followups/task1_sync.sh both 2>&1 | tail -2
  sleep 600
done
