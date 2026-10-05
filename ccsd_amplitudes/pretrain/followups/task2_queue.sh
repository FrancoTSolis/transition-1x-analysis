#!/bin/bash
# Run task-2 optimization jobs sequentially: one "args for task2_opt" per line of $1 (lines starting with # skipped).
# A job is skipped when its result JSON exists.  GPU: scai2 GPU 3 only (CUDA_VISIBLE_DEVICES=3).
# usage: setsid nohup bash pretrain/followups/task2_queue.sh <jobs.txt> > <log> 2>&1 < /dev/null &
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
PY=pretrain/.train_venv/bin/python3
export CUDA_VISIBLE_DEVICES=3
RES=pretrain/opt_true/results/followups/task2_direct_opt
while IFS= read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  set -- $line
  name= start= obj=lucj opt=nomad tag=
  while [[ $# -gt 0 ]]; do
    case "$1" in --name) name=$2;; --start) start=$2;; --objective) obj=$2;; --optimizer) opt=$2;; --tag) tag=$2;; esac
    shift
  done
  stem="${name}__${start}__${opt}"; [[ -n "$tag" ]] && stem="${stem}__${tag}"
  if [[ -f "$RES/$obj/$stem.json" ]]; then echo "skip $obj/$stem"; continue; fi
  echo "=== $(date '+%F %T') start $obj/$stem"
  $PY -m pretrain.followups.task2_opt $line 2>&1 | grep -v -E "^[0-9]+ +-?[0-9.]+ *$|Best feasible|^Warning" | tail -4
  echo "=== $(date '+%F %T') end $obj/$stem"
done < "$1"
echo "queue done $(date '+%F %T')"
