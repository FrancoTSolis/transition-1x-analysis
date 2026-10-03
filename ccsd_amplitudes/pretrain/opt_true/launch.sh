#!/bin/bash
# Launch one pretrain.opt_true.train run detached on a given host/GPU (shared filesystem).
# Usage: launch.sh <host|local> <gpu> <run_name> [train.py args...]
#   e.g. launch.sh scai4 3 full_long_big --frame canonical --epochs 300 ...
# Logs to runs_ot/<run_name>.log; outputs to runs_ot/<run_name>/.
set -euo pipefail
HOST=$1; GPU=$2; NAME=$3; shift 3
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
ARGS=$(printf ' %q' "$@")
CMD="cd $ROOT; CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=3 MKL_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 setsid nohup pretrain/.train_venv/bin/python3 -u -m pretrain.opt_true.train --out runs_ot/$NAME$ARGS > runs_ot/$NAME.log 2>&1 < /dev/null &"
if [ "$HOST" = "local" ] || [ "$HOST" = "$(hostname -s)" ]; then
  bash -c "$CMD"
else
  ssh -o BatchMode=yes -o ConnectTimeout=10 "fts@$HOST.cs.ucla.edu" "$CMD" < /dev/null
fi
echo "launched $NAME on $HOST:GPU$GPU"
