#!/bin/bash
# Task 5: start one detached dm-engine benchmark process (pretrain/followups/task5_dm_bench.py) on THIS host.
#   bash pretrain/followups/task5_launch.sh <gpu> <label> <mem_cap_gb> <jobs> [extra bench args...]
# Run it on the GPU host (e.g. through ssh scai7).  4 CPU threads per process (block2 <H> + host-side numpy).
# Refuses to start if the user-GPU-limit guard (<= 2 GPUs per machine on scai3-7) says no.
set -euo pipefail
GPU=$1; LABEL=$2; CAP=$3; JOBS=$4; shift 4
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
PY=$ROOT/pretrain/.train_venv/bin/python3
LOG=$ROOT/rl_runs/followups/task5_dm_scaling
OUT=$ROOT/pretrain/opt_true/results/followups/task5_dm_scaling/bench.jsonl
cd $ROOT
$PY -m pretrain.rl.reward_queue guard --gpu $GPU || exit 3
export CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 RAYON_NUM_THREADS=4 \
       NUMBA_NUM_THREADS=4
setsid nohup $PY -u pretrain/followups/task5_dm_bench.py --params rl_runs/followups/task5_dm_scaling/params_rl4n29f.pkl \
    --bench-params C3H5N3_rxn2003_P:net:rl_runs/followups/task5_dm_scaling/bench_params_C3H5N3_rxn2003_P.npz \
    --jobs "$JOBS" --b2-threads 4 --mem-cap-gb $CAP --label $LABEL --out $OUT "$@" \
    > $LOG/bench_$LABEL.log 2>&1 < /dev/null &
echo "$LABEL pid $! on $(hostname -s) GPU $GPU"
