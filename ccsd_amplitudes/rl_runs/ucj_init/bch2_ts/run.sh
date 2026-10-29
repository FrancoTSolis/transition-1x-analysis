#!/bin/bash
# BCH2 (Avijit's optimizer, unmodified) on two norb-15 TS-set molecules from 5 seeds; 4 jobs at a time, 4 threads each
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
jobs=()
for name in C2H3N_rxn2858_TS C2H3N_rxn2858_P; do
  for seed in default rl4n29f label truncated random; do
    jobs+=("$name $seed")
  done
done
printf '%s\n' "${jobs[@]}" | xargs -P 4 -L 1 bash -c 'pretrain/.train_venv/bin/python -m ucj_init.bch2_ts --name $0 --seed $1 --max-cycle 25 --out pretrain/opt_true/results/ucj_init/bch2_ts/$0__$1.json > rl_runs/ucj_init/bch2_ts/$0__$1.log 2>&1'
