#!/bin/bash
# Task 7: one more seed of the norb 15-18 exact-reward GRPO run (recipe of rl_runs/grpo_large_rl4/args.json).
# Usage: bash pretrain/followups/task7_launch_seed.sh <seed>     (from ccsd_amplitudes/, on scai2)
# Only differences from the seed-0 run: --seed, --queue-root rl_queue/exact_shared, --queue-timeout 21600
# (the shared queue also serves other agents; a 1-h timeout would turn a long wait into NaN rewards),
# --driver-threads 4 and explicit thread caps (CPU budget).
set -euo pipefail
SEED=$1
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
OUT=rl_runs/followups/task7_rl_seeds/seed${SEED}
if [ -e "$OUT/log.jsonl" ]; then echo "$OUT/log.jsonl exists; refusing to append" >&2; exit 1; fi
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 RAYON_NUM_THREADS=4 NUMBA_NUM_THREADS=4
setsid nohup pretrain/.train_venv/bin/python3 -u -m pretrain.rl.grpo_slot \
    --init runs_ot/slotall_T4/best.pt --prefix-T 3 --init-policy rl_runs/grpo_slot_T4p/policy_best.pt \
    --train-names pretrain/rl/large_train.txt --val-names pretrain/rl/large_val.txt \
    --reward queue --queue-root rl_queue/exact_shared --queue-timeout 21600 \
    --group-size 8 --batch-mols 4 --sigma-k 0.01 --sigma-z 0.005 --antithetic --ppo-epochs 1 --lr 2e-5 \
    --steps 50 --eval-every 5 --device cpu --driver-threads 4 --seed ${SEED} --out ${OUT} \
    > rl_runs/followups/task7_rl_seeds/seed${SEED}.log 2>&1 < /dev/null &
echo "seed ${SEED}: pid $!"
