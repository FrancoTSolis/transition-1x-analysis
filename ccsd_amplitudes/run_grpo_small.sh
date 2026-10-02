#!/bin/bash
# GRPO fine-tuning PoC: 30 smallest molecules (norb<=16), exact statevector reward,
# policy initialised from the amortized compressed-DF model (backbone + residual heads).
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=${GPU:-1}
OUT=${OUT:-rl_runs/small_exact_v1}
mkdir -p $OUT
exec taskset -c ${CORES:-28-47} pretrain/.train_venv/bin/python3 -u -m pretrain.rl.grpo \
    --names-file gauge_study/names_small_norb_le16.txt --val-frac 0.2 \
    --connectivity square --reward exact \
    --init-backbone ${INIT:-checkpoints_comp_recon_n29/best.pt} \
    --group-size ${G:-8} --batch-mols ${B:-4} --sigma ${SIGMA:-0.02} \
    --clip-eps 0.2 --kl-coef 0.0 --ref-update 5 --anchor-coef 0.0 \
    --lr ${LR:-2e-5} --steps ${STEPS:-60} --eval-every 10 \
    --n-workers ${NW:-20} --worker-threads 1 --out-dir $OUT
