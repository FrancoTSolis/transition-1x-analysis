#!/bin/bash
# Resume (5 Oct, successor agent): score each new objective-A NOMAD result with the Task-1 protocol as soon as its JSON
# appears (3 SCI threads, scai2 GPU 3); exits after the reordered NOMAD queue (chain_nomad2.sh) is done and all are scored,
# or after 3 attempts in a row that left work undone (first launch hit a CUDA OOM next to a 6 GB dimscan worker).
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
D=rl_runs/followups/task2_direct_opt
R=pretrain/opt_true/results/followups/task2_direct_opt
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
fails=0; last=
while true; do
  done_queue=0; grep -q "nomad2 queue done" $D/queue_A_nomad2.log && done_queue=1
  todo=$(pretrain/.train_venv/bin/python3 - <<'PY'
import json, glob
R = "pretrain/opt_true/results/followups/task2_direct_opt"
have = {(json.loads(l)["name"], json.loads(l)["cand"]) for l in open(f"{R}/score_task1_protocol.jsonl")}
names = set()
for f in glob.glob(f"{R}/lucj/*__nomad.json"):
    r = json.load(open(f))
    if (r["name"], f"lucj:{r['start']}:nomad") not in have:
        names.add(r["name"])
print(" ".join(sorted(names)))
PY
)
  if [ -n "$todo" ] && [ "$todo" = "$last" ]; then fails=$((fails + 1)); else fails=0; fi
  if [ $fails -ge 3 ]; then echo "=== $(date '+%F %T') giving up on: $todo"; break; fi
  last=$todo
  if [ -n "$todo" ]; then
    echo "=== $(date '+%F %T') scoring $todo"
    pretrain/.train_venv/bin/python3 -m pretrain.followups.task2_score --only lucj --threads 3 --max-mem-gb 2.0 \
        --names $todo < /dev/null 2>&1 | grep --line-buffered -v -i warn
  elif [ $done_queue = 1 ]; then
    echo "=== $(date '+%F %T') follow done"; break
  else
    sleep 60
  fi
done
