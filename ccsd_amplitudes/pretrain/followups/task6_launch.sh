#!/bin/bash
# Task 6: n29 GRPO with dm-engine chi-128 TN rewards on scai7 GPU 4 (4 workers x 4 threads), driver on scai2 CPU.
#   bash pretrain/followups/task6_launch.sh workers      # (run on scai7) 4 TN chi-128 dm workers on GPU 4
#   bash pretrain/followups/task6_launch.sh driver S B   # (run on scai2) GRPO driver: S steps, B molecules per step
#   bash pretrain/followups/task6_launch.sh v1workers    # (run on scai7) 4 v1 chi-256 evaluation workers on GPU 4
# Recipe = rl_runs/grpo_n29_tn/args.json (Expanse run) except: the reward (dm chi 128 through the queue instead of
# v1 chi 64 in-process), --batch-mols / --steps (chosen from the measured throughput), --queue-timeout 21600,
# --driver-threads 4.
set -euo pipefail
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
cd $ROOT
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 RAYON_NUM_THREADS=4 NUMBA_NUM_THREADS=4
case $1 in
  workers)
    [ "$(hostname -s)" = scai7 ] || { echo "run on scai7"; exit 1; }
    for s in a b c d; do
      THREADS=4 WORKER_ARGS="--kind tn --chi 128 --tn-impl current --block2-threads 4" \
        bash pretrain/rl/launch_workers.sh scai7 4 $ROOT/rl_queue/task6_tn128 44 complex64 "" 1 $s
    done ;;
  v1workers)
    [ "$(hostname -s)" = scai7 ] || { echo "run on scai7"; exit 1; }
    for s in a b c d; do
      THREADS=4 WORKER_ARGS="--kind tn --chi 256 --tn-impl v1 --block2-threads 4" \
        bash pretrain/rl/launch_workers.sh scai7 4 $ROOT/rl_queue/task6_v1_256 29 complex64 "" 1 $s
    done ;;
  driver)
    [ "$(hostname -s)" = scai2 ] || { echo "run on scai2"; exit 1; }
    STEPS=$2; BM=$3
    OUT=rl_runs/followups/task6_n29_rl_dm/run
    if [ -e "$OUT/log.jsonl" ]; then echo "$OUT/log.jsonl exists; refusing to append" >&2; exit 1; fi
    mkdir -p $OUT
    setsid nohup pretrain/.train_venv/bin/python3 -u -m pretrain.rl.grpo_slot \
        --init runs_ot/slotall_T4/best.pt --prefix-T 3 --init-policy rl_runs/grpo_large_rl4/policy_snap_for_n29.pt \
        --train-names pretrain/rl/n29_train40.txt --val-names pretrain/rl/n29_val20.txt \
        --reward queue --queue-root rl_queue/task6_tn128 --queue-kind tn --queue-timeout 21600 \
        --group-size 8 --batch-mols $BM --sigma-k 0.01 --sigma-z 0.005 --antithetic --ppo-epochs 1 --lr 1e-5 \
        --steps $STEPS --eval-every 5 --device cpu --driver-threads 4 --seed 0 --out $OUT \
        > rl_runs/followups/task6_n29_rl_dm/run.log 2>&1 < /dev/null &
    echo "driver pid $!" ;;
  *) echo "usage: $0 workers|v1workers|driver S B"; exit 1 ;;
esac
