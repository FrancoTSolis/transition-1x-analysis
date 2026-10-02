#!/bin/bash
# Supervised regression of canonical-gauge compressed labels (dkappa, dZ residuals), n29 group.
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=${GPU:-0}
CFG=${CFG:-square_reg0.005}
exec pretrain/.train_venv/bin/python3 -m pretrain.train --loss-mode compressed_sup \
    --comp-dir rhf_targets_compressed/$CFG --connectivity square \
    --names-file gauge_study/names_n29_16_13.txt \
    --epochs 40 --batch-size 12 --lr 5e-4 --warmup-steps 200 \
    --embed-dim 192 --num-layers 6 --num-heads 8 --log-interval 20 --num-workers 3 \
    --checkpoint-dir checkpoints_comp_sup_n29_$CFG
