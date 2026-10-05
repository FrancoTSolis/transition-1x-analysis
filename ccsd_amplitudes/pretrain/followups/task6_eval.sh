#!/bin/bash
# Task 6, after the RL run: v1 chi-256 energies (the docs' evaluation setting) of the run's best and last policies on
# n29 val20 / test24, plus a 4-energy control (stored rl4n29f v1 chi-256 energies recomputed on this GPU).
# Run on scai2 once the v1 workers serve rl_queue/task6_v1_256 (task6_launch.sh v1workers on scai7).
#   bash pretrain/followups/task6_eval.sh
set -euo pipefail
ROOT=/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
cd $ROOT
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 RAYON_NUM_THREADS=4 NUMBA_NUM_THREADS=4
PY=pretrain/.train_venv/bin/python3
RUN=rl_runs/followups/task6_n29_rl_dm/run
TD=rl_runs/followups/task6_n29_rl_dm
OUTD=pretrain/opt_true/results/followups/task6_n29_rl_dm
Q=rl_queue/task6_v1_256
SB=$($PY -c "import torch; print(torch.load('$RUN/policy_best.pt', map_location='cpu', weights_only=False)['step'])")
SL=$($PY -c "import torch; print(torch.load('$RUN/policy_last.pt', map_location='cpu', weights_only=False)['step'])")
echo "best step $SB, last step $SL"
POL="--policy $RUN/policy_best.pt --tag rl4n29dm"
[ "$SB" != "$SL" ] && POL="$POL --policy $RUN/policy_last.pt --tag rl4n29dmf"
for S in val20:n29val test24:n29test; do
  L=${S%%:*}; T=${S##*:}
  $PY -m pretrain.rl.policy_dump --names-file pretrain/rl/n29_$L.txt $POL --device cpu --out $TD/tasks_${T}_rl4n29dm.pkl
done
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
$PY -u -m pretrain.rl.queue_eval --kind tn --dump runs_ot/energy_tasks/n29val_rl4n29.pkl \
    --only-names $(cat $TD/control_val4.txt) --skip-keys rl4n29 --queue-root $Q \
    --out $OUTD/energy_n29val_control_v1chi256.json > $TD/eval_control.log 2>&1 &
$PY -u -m pretrain.rl.queue_eval --kind tn --dump $TD/tasks_n29val_rl4n29dm.pkl --queue-root $Q \
    --out $OUTD/energy_n29val_rl4n29dm_v1chi256.json > $TD/eval_val.log 2>&1 &
$PY -u -m pretrain.rl.queue_eval --kind tn --dump $TD/tasks_n29test_rl4n29dm.pkl --queue-root $Q \
    --out $OUTD/energy_n29test_rl4n29dm_v1chi256.json > $TD/eval_test.log 2>&1 &
wait
echo "evaluations done $(date)"
