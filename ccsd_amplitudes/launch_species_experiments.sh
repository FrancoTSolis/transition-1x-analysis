#!/bin/bash
set -euo pipefail

VENV=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes/pretrain/.train_venv/bin/python3
WORKDIR=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
COMMON_ARGS="--jobs-dir jobs --epochs 200 --batch-size 4 --lr 1e-4 --warmup-steps 200 --log-interval 5 --embed-dim 128 --num-layers 4 --num-heads 8 --no-amp --dropout 0.0"

cd "$WORKDIR"
mkdir -p checkpoints_ts checkpoints_rp

echo "Launching TS-only on GPU 3..."
CUDA_VISIBLE_DEVICES=3 nohup $VENV -m pretrain.train $COMMON_ARGS \
    --species-filter TS \
    --checkpoint-dir checkpoints_ts \
    > pretrain_ts_stdout.log 2>&1 &
echo "  PID=$!"

echo "Launching RP-only on GPU 4..."
CUDA_VISIBLE_DEVICES=4 nohup $VENV -m pretrain.train $COMMON_ARGS \
    --species-filter RP \
    --checkpoint-dir checkpoints_rp \
    > pretrain_rp_stdout.log 2>&1 &
echo "  PID=$!"

echo "Species experiments launched."
