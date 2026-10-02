#!/bin/bash
# Canonical compressed-DF labels for the n29 (16,13) group on scai1: 60 pinned single-core shards.
# Configs run sequentially inside each shard so the first config finishes first.
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
mkdir -p gen_logs/compressed
N=60
for i in $(seq 0 $((N-1))); do
  core=$((i + 10))   # leave cores 0-9 for other users
  nohup taskset -c $core .lucj_venv/bin/python3 -u generate_compressed_targets.py \
      --names-file gauge_study/names_n29_16_13.txt --out-dir rhf_targets_compressed \
      --configs square_reg0.005 square_reg0 all-to-all_reg0.005 --maxiter 500 \
      --shard $i --n-shards $N > gen_logs/compressed/shard_$i.log 2>&1 &
done
wait
