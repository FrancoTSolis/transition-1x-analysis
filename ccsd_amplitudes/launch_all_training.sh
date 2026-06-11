#!/bin/bash
set -euo pipefail

VENV=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes/pretrain/.train_venv/bin/python3
WORKDIR=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
COMMON_ARGS="--jobs-dir jobs --epochs 200 --batch-size 4 --lr 1e-4 --warmup-steps 200 --log-interval 5 --embed-dim 128 --num-layers 4 --num-heads 8 --no-amp --dropout 0.0"

cd "$WORKDIR"
rm -rf checkpoints checkpoints_hf checkpoints_mo
mkdir -p checkpoints checkpoints_hf checkpoints_mo

echo "Launching positional on GPU 0..."
CUDA_VISIBLE_DEVICES=0 nohup $VENV -m pretrain.train $COMMON_ARGS \
    --checkpoint-dir checkpoints \
    > pretrain_pos_stdout.log 2>&1 &
echo "  PID=$!"

echo "Launching HF energies on GPU 1..."
CUDA_VISIBLE_DEVICES=1 nohup $VENV -m pretrain.train $COMMON_ARGS \
    --use-hf-energies \
    --checkpoint-dir checkpoints_hf \
    > pretrain_hf_stdout.log 2>&1 &
echo "  PID=$!"

echo "Launching MO coefficients on GPU 2..."
CUDA_VISIBLE_DEVICES=2 nohup $VENV -m pretrain.train $COMMON_ARGS \
    --use-mo-coeffs \
    --checkpoint-dir checkpoints_mo \
    > pretrain_mo_stdout.log 2>&1 &
echo "  PID=$!"

echo "All 3 launched."
