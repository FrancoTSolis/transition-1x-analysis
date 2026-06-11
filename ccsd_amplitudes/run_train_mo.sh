#!/bin/bash
set -euo pipefail

VENV=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes/pretrain/.train_venv/bin/python3
WORKDIR=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes

cd "$WORKDIR"
mkdir -p checkpoints_mo

export CUDA_VISIBLE_DEVICES=2

exec "$VENV" -m pretrain.train \
    --jobs-dir jobs \
    --epochs 200 \
    --batch-size 4 \
    --lr 1e-4 \
    --warmup-steps 200 \
    --log-interval 5 \
    --embed-dim 128 \
    --num-layers 4 \
    --num-heads 8 \
    --no-amp \
    --dropout 0.0 \
    --use-mo-coeffs \
    --checkpoint-dir checkpoints_mo
