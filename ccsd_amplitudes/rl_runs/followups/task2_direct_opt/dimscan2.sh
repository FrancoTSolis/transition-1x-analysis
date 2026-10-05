#!/bin/bash
# Resume (5 Oct, successor agent): fixed-size subspace scan of two norb-16 molecules, label / one-call starts and their
# objective-A SPSA optima (does the LUCJ-energy optimization pick better configurations at fixed size?).  3 threads, GPU 3.
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=1 NUMBA_NUM_THREADS=1
echo "=== $(date '+%F %T') dimscan (norb 16)"
pretrain/.train_venv/bin/python3 -u -m pretrain.followups.task2_dimscan --names C2H4O_rxn0724_TS C3H4_rxn2391_TS --only lucj \
    --cands start:label start:rl4n29f lucj:label:spsa lucj:rl4n29f:spsa:aLabel --threads 3 --max-mem-gb 2.5 \
    < /dev/null 2>&1 | grep --line-buffered -v -i warn
echo "=== $(date '+%F %T') dimscan done"
