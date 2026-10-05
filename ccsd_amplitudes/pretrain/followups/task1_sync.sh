#!/bin/bash
# Task 1: push new sample files to Expanse and pull back the diagonalization results of the Expanse workers.
# Expanse only ever writes the units it was given (disjoint from the local workers'), so the pull cannot clobber
# local results; claim files and partial *.tmp files are excluded.
#   bash pretrain/followups/task1_sync.sh [push|pull|both]
set -euo pipefail
CA=$(cd "$(dirname "$0")/../.." && pwd)
cd "$CA"
R=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
W=runs_ot/energy_tasks/task1_sqd_vs_lin
EX="python3 expanse/expanse.py"
mode=${1:-both}
if [[ $mode == push || $mode == both ]]; then
  $EX put pretrain/followups/task1_sqd_lin.py $R/pretrain/followups/task1_sqd_lin.py
  $EX put $W/samples/ $R/$W/samples/ --exclude='*.part'
fi
if [[ $mode == pull || $mode == both ]]; then
  $EX get $R/$W/diag/ $W/diag/ --exclude=claims/ --exclude='*.tmp'
  $EX get $R/rl_runs/followups/task1_sqd_vs_lin/ rl_runs/followups/task1_sqd_vs_lin/expanse_logs/
fi
echo "local diag units: $(ls $W/diag/*.json 2>/dev/null | wc -l)"
