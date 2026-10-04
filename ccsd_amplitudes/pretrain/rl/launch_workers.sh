#!/bin/bash
# Start one reward-queue worker on <host> GPU <gpu> (shared filesystem; logs in <queue_root>/logs/).
#   launch_workers.sh <host|local> <gpu> <queue_root> <max_norb> [dtype] [max_mem_gb] [max_items] [id_suffix]
#   extra worker flags via the WORKER_ARGS environment variable, e.g. WORKER_ARGS="--kind tn --chi 256";
#   CPU threads per worker via THREADS (default 4; a TN worker uses about THREADS cores while it runs)
# The worker runs the user-GPU-limit guard itself (at most 2 GPUs per machine on scai3-scai7, counting the user's
# other jobs) at start-up and every minute while running; this script also refuses up front.
set -euo pipefail
HOST=$1; GPU=$2; Q=$3; MAXN=$4; DT=${5:-complex64}; MEM=${6:-}; ITEMS=${7:-1}; SUF=${8:-}
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
PY=$ROOT/pretrain/.train_venv/bin/python3
WID="${HOST}-gpu${GPU}${SUF:+-$SUF}"
TH=${THREADS:-4}
MEMARG=""; [ -n "$MEM" ] && MEMARG="--max-mem-gb $MEM"
CMD="cd $ROOT; mkdir -p $Q/logs; $PY -m pretrain.rl.reward_queue guard --gpu $GPU || exit 3; CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=$TH MKL_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=$TH setsid nohup $PY -u -m pretrain.rl.reward_queue worker --root $Q --worker-id $WID --max-norb $MAXN --dtype $DT --max-items $ITEMS $MEMARG ${WORKER_ARGS:-} > $Q/logs/$WID.log 2>&1 < /dev/null &"
if [ "$HOST" = "local" ] || [ "$HOST" = "$(hostname -s)" ]; then
  bash -c "$CMD"
else
  timeout 60 ssh -n -o BatchMode=yes -o ConnectTimeout=10 "fts@$HOST.cs.ucla.edu" "$CMD"
fi
echo "launched worker $WID (max_norb $MAXN, $DT) -> $Q"
