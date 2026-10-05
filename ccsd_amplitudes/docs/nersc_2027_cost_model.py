#!/usr/bin/env python3
"""ML-LUCJ GPU budget model for the 2027 NERSC (Perlmutter) renewal.  Transparent, runnable, stdlib only.
Version 2 (4 Oct 2026): revised after the arithmetic / provenance / reviewer checks (see CHANGES below).

Unit of the request: Perlmutter GPU node-hour (1 node for 1 h = 4 x NVIDIA A100; charge factor 1, regular QOS).
Single-GPU jobs in the "shared" QOS are charged by GPU fraction, so node-hours = A100-hours / 4 throughout.

Every input is a named constant.  Scenario inputs are T(low, central, high); "low" is always the value that makes
the budget SMALLER (so for a speed-up factor low > high).  Reported scenarios:
  CENTRAL   every input at its central value (the planning number)
  LOW/HIGH  Monte-Carlo P10 / P90 per line item and for the total (4,000 draws, each input log-triangular between
            its low and high ends with mode = central); percentiles are per line and do not add up
  corners   every input at its cheap / expensive end at once (printed for transparency; not a planning range)
Source tags:
  [measured: <file>]   our own logs; paths relative to REPO
  [derived: ...]       arithmetic or extrapolation from measured values
  [paper: <ref>]       Lin et al. arXiv:2511.22476 (local copy wan-hsuan-lucj/) or another publication / web page
  [inferred: ...]      a reading of a source that the source does not state
  [assumption: <why>]  a planning choice
Usage:  python3 cost_model.py         unit costs, line items (central, P10/P50/P90, corners), sensitivity
        python3 cost_model.py --md    also print the markdown tables used in the note

CHANGES v1 -> v2
  - SQD_DIM central 1e7 -> 4e6 (range 1e6-1.6e7): Lin et al.'s simulated subspaces never reached the 4,000^2 cap;
    the v1 'measured' justification came from unfaithful MPS samples.
  - SBD_NLINK_EXP central 1 -> 2 (range 1-2.5): with exponent 1 the GPU solver would be 56-126x faster per node than
    the CPU path, but the published like-for-like speedup is ~17-29x per Perlmutter node; exponent 2 reconciles our
    measured CPU cost at norb 29, the published GPU anchors (nlink ~100-110) and that speedup (implied ratio 27x).
  - A_STATE central 2.5 -> 2.6 (= the chi-128 fit); B_CHI_STATE high end 2.5 -> 3.0; LAB_TO_A100 high end 1.0 -> 0.8.
  - New line: training labels for the orbital-ordering model (project outline, part 1).
  - New outputs: implied GPU/CPU node ratio, request percentiles, per-molecule TN-optimization range, reduced-
    max_dim CPU fallback, storage, formulas with explicit unit conversions.
"""
import math
import random
import statistics
import sys

REPO = "/xuanwu-tank/east/fts/projects/transition-1x-analysis/"
A100_PER_NODE = 4            # [paper: docs.nersc.gov/systems/perlmutter/architecture] 4x A100 per GPU node
CORES_PER_GPU_NODE = 64      # [same] one AMD EPYC 7763 per GPU node -> 16 cores per A100
CORES_PER_CPU_NODE = 128     # [same] two EPYC 7763 per CPU node
SCEN = ("low", "central", "high")


def T(lo, c, hi):
    return (lo, c, hi)


# =================================================================================================================
# 0. Dataset and splits  [measured: ccsd_amplitudes/rhf_dataset/_index.json]  30,205 molecules (norb -> count)
# =================================================================================================================
NORB_COUNTS = {15: 6, 16: 24, 17: 21, 18: 33, 19: 12, 20: 47, 21: 48, 22: 144, 23: 276, 24: 375, 25: 336, 26: 624,
               27: 572, 28: 1602, 29: 2055, 30: 2766, 31: 2436, 32: 3493, 33: 1690, 34: 3112, 35: 2224, 36: 3204,
               37: 1578, 38: 1997, 39: 335, 40: 947, 41: 57, 42: 182, 44: 9}
# median number of occupied (alpha) orbitals per norb [measured: same file]; nlink = nocc * (nvirt + 1)
NOCC_MED = {15: 8, 16: 9, 17: 10, 18: 10, 19: 11, 20: 11, 21: 12, 22: 12, 23: 13, 24: 13, 25: 14, 26: 15, 27: 16,
            28: 16, 29: 16, 30: 17, 31: 17, 32: 17, 33: 18, 34: 19, 35: 19, 36: 19, 37: 20, 38: 20, 39: 21, 40: 21,
            41: 21, 42: 21, 44: 22}
# Proposed formula-disjoint split: official Transition1x split + norb 39-44 holdout (C5H11NO, C6H14O, C7H16)
# [measured: data/transition1x.h5 split groups x rhf_dataset/_index.json formulas]
TEST = {17: 12, 25: 33, 29: 486, 31: 15, 32: 24, 33: 291, 39: 242, 42: 62, 44: 9}   # 1,174 molecules, 11 formulas
VAL = {21: 9, 22: 42, 24: 12, 26: 6, 30: 288, 37: 303, 38: 15}                        # 675 molecules, 8 formulas
TRAIN = {n: NORB_COUNTS[n] - TEST.get(n, 0) - VAL.get(n, 0) for n in NORB_COUNTS}     # 28,356 molecules
BUCKETS = {"B1": (15, 19), "B2": (20, 29), "B3": (30, 37), "B4": (38, 44)}
# Central SQD/QSD test subset: 300 test molecules (~100 reactions x R/P/TS) over all 11 test formulas
# [assumption: balanced reporting; C3H5NO2 + C5H5NO are 66% of the test set, so they are subsampled]
SQD_SUBSET = {17: 12, 25: 33, 29: 60, 31: 15, 32: 24, 33: 60, 39: 60, 42: 27, 44: 9}
N_ALL = sum(NORB_COUNTS.values())

# =================================================================================================================
# 1. Measured unit costs (fixed)
# =================================================================================================================
# Exact LUCJ energy on one A100 (pretrain/rl/gpu_energy.py, complex64), A100-seconds per energy.
# [measured: rl_queue/logs/bench_big_scai3.log, A100-40GB: norb17 (10,7) 3.0 s; norb18 (11,7) 9.2, (10,8) 17.9,
#  (9,9) 22.7 s -> weighted by the dataset shapes 18.0 s.  derived: norb15-16 and norb17 (9,8) are TITAN Xp times
#  (pretrain/rl/tests/results/bench_gpu_energy_titanxp.log: 0.66 s, 2.2-2.9 s, 11.5 s) / 2.44 (A100/TITAN Xp ratio
#  measured at norb 17-18).  norb19 (11,8): H200 shared with another user's job, 43-57 s at 55 GB
#  (pretrain/opt_true/results/energy_n19_all.json) x H200/A100 ratio 1.6-2.2 -> ~100 s; needs an 80 GB A100 node]
EXACT_A100_S = {15: 0.27, 16: 1.0, 17: 3.8, 18: 18.0, 19: 100.0}

# Tensor-network (MPS) LUCJ energy, density-matrix engine pretrain/rl/tn_energy.py (method="dm", the RL reward)
# [measured: pretrain/rl/tests/results/tn_v2/*.jsonl, medians at norb 29 (58 qubits); single worker on a TITAN Xp
#  on a heavily loaded host (load average 60-120 on 48 cores); block2 <H> on 1 CPU thread]
T_STATE_N29_TXP = {64: 65.7, 128: 151.1, 256: 536.3}       # s, MPS build (GPU + host-side fp64 density matrices)
T_EXPECT_N29_1THR = {64: 57.0, 128: 123.1, 256: 289.6}     # s, block2 <H> (CPU)
# chi 512: ONE run (C3H5N3_rxn2003_P, n29_dm.jsonl, host load 60): state 774 s vs 309 s, <H> 628 s vs 205 s at
# chi 256 for the same molecule (load 85) -> scale the chi-256 medians by these same-molecule ratios (x2.51, x3.06)
T_STATE_N29_TXP[512] = 536.3 * 774.5 / 309.0
T_EXPECT_N29_1THR[512] = 289.6 * 628.0 / 205.0
T_MPO_N29 = 40.0             # s, block2 MPO build, once per molecule [measured: 33-49 s, tn_v2 t_mpo_build]
B2_4THR_SPEEDUP = 126.1 / 36.3   # [measured: tn_v2/expect_threads_n29.log, chi 128: 1 thread 126 s, 4 threads 36 s]
WORKERS_PER_GPU = 4          # [assumption: CPU-core limited, 16 cores per A100 / 4 block2 threads; GPU memory per
                             #  worker 1.3 / 3.7-4.0 GB torch peak at chi 128 / 256 (2-3x that in reality); chi-256
                             #  evaluations at norb >= 37 go to the 80 GB nodes (same charge); 4 workers per A100
                             #  assume the MPS build is about half GPU work (not yet measured on an A100)]
A_MPO = 4.0                  # [assumption: number of Hamiltonian terms ~ n^4]
EXACT_MAX_NORB_RL = 18       # [assumption: exact GPU rewards up to norb 18 (fits 40 GB); norb 19 needs an 80 GB
                             #  node and costs ~100 A100-s vs ~10 for a chi-128 TN energy, so use TN there]
GROUP = 8                    # GRPO group size [measured: rl_runs/grpo_n29_tn/args.json]
DEV_RUN_ENERGIES = 30 * 12 * 8 + 220   # one n29-scale development run [measured: rl_runs/grpo_n29_tn: 30 steps x
                                       #  12 molecules x 8 samples + 220 evaluation energies]
V1_CHI64_CPU_ENERGIES_PER_NODE_H = 96 * 3600 / 262.0   # CPU fallback for rewards [measured: rl_runs/grpo_n29_tn/
                                       #  log.jsonl, 96 single-thread workers, 262 s per 96-energy step, Expanse node]

# SQD protocol of Lin et al. [paper: main.tex:322 (1e5 samples, 10 batches x 4,000, QSCI), 325-332 (MPS chi 50,
#  1e4 samples per objective, NOMAD 500 evaluations), 335-336 (10 runs x 1e6 shots, 10 configuration-recovery
#  iterations); max_dim 4,000 per spin -> <= 1.6e7 determinants: scripts/paper/n2_cc-pvdz_10e26o/tn_results.py:49]
N_BATCH = 10
TN_SQD_ITERS = 1             # simulation / TN evaluation: one QSCI pass per batch
QPU_CR_ITERS = 10            # hardware: 10 configuration-recovery iterations
QPU_D = 1.6e7                # hardware subspaces assumed at the max_dim cap (noisy samples + recovery fill it)
QPU_NORB = 29                # QPU molecules: test formula C3H5NO2, 58 data qubits (+ ancillas on heavy-hex)
NLINK_N2 = 5 * (21 + 1)      # [inferred: Walkup et al. do not state the basis of their N2 input; assumed the
                             #  (10e,26o) cc-pVDZ N2 of IBM's SQD papers -> nocc 5, nvirt 21, nlink 110.  The H2O
                             #  anchor (cc-pVDZ, 24 orbitals) has nlink 100.  Replace with our own SBD benchmark]
BASE_EVALS, BASE_SHOTS, BASE_BATCHES, BASE_D, BASE_CHI, BASE_NORB = 500, 1e4, 3, 1e6, 128, 29
# per-instance TN-optimization baseline: Lin et al.'s NOMAD loop (500 evaluations, 1e4 MPS samples each) but with
# reduced SQD settings (3 batches, max_dim 1000) [assumption: full settings cost hundreds of node-h per molecule]
# Published GPU-vs-CPU speedup of the SBD solver, for the consistency check of the two diagonalization paths:
# [paper: Walkup et al. arXiv:2601.16169 Table II: H2O 6.28e8 dets on 32 Frontier nodes, GPU 58.3 s vs CPU
#  3,652 s (cached) / 5,568 s -> 63x / 95x per Frontier node; Table IV: 8 A100 309 s vs 16 MI250X GCDs 190 s
#  -> 1 A100 = 1.23 GCD].  Frontier node = 8 GCDs + one 64-core CPU (56 cores usable); Perlmutter GPU node =
#  4 A100 = 4.9 GCDs; Perlmutter CPU node = 128 cores.
WALKUP_SPEEDUP_PER_FRONTIER_NODE = (63.0, 95.0)
A100_IN_GCD = 3040.0 / 2472.0
FRONTIER_CPU_CORES = (56, 64)
# Orbital-ordering model (project outline NN_LUCJ_full_outline, part 1 "site ordering"): one training label =
# molecule + hardware topology + layout + LUCJ energy from the MPS judge (ordering_ml_diagrams.tex); the O3 pilot
# used 15 orderings x 2 topologies (outline slide 3)
ORD_TOPOLOGIES = 2           # heavy-hex + square

# =================================================================================================================
# 2. Scenario inputs  T(low, central, high)   ("low" = cheaper)
# =================================================================================================================
PARAMS = dict(
    # ---- Pretraining (label-free compressed-DF objective, slot transformer) ----
    EPOCH_LAB_GPU_H=T(0.18, 0.30, 0.45),   # [measured: runs_ot/slotall_{T1,T8,T4,distill_T1}/log.jsonl] GPU-h per
                                           #  epoch over 30,034 molecules, 0.92M-param model: 0.18 / 0.23 / 0.43 /
                                           #  0.45 (4 runs); GPU type not logged (lab: TITAN Xp, V100, A100, A6000,
                                           #  L40S); the 3.0M model with 2x batches: 0.08 (runs_ot/slotall_big_T1)
    LAB_TO_A100=T(2.5, 1.5, 0.8),          # [measured 2.4-2.5x A100/TITAN Xp on the exact kernel; assumption:
                                           #  training is overhead-bound; 0.8 = the lab card was faster than an A100]
    N_PROD=T(6, 10, 12),                   # [assumption] final runs: one-shot + recycled, square + heavy-hex, CCSD/MP2 inputs, seeds
    EPOCHS_PROD=T(150, 300, 400),          # [assumption] long schedules: our fitting studies needed hundreds of epochs
    MODEL_MULT_PROD=T(4.0, 8.0, 12.0),     # [assumption] production models up to ~40M params (pair-token transformers), sub-linear
    N_SWEEP=T(100, 160, 240),               # [assumption] architectures (slot / pair transformer, sparse attention,
                                           #  recycles, output parameterizations, amplitude compression) x hyper-
                                           #  parameters x 2 seeds; the Oct-2 campaign alone was 51 runs in a day
    EPOCHS_SWEEP=T(20, 40, 60),            # [assumption] full-data-epoch equivalents per sweep config
    MODEL_MULT_SWEEP=T(2.5, 4.0, 6.0),     # [assumption] sweeps reach d256-d512 (5-42M params)
    # ---- TN energy engine scaling ----
    A100_STATE_SPEEDUP=T(2.5, 2.0, 1.5),   # [assumption: A100 vs loaded TITAN Xp; 2.4-2.5x measured on the exact
                                           #  kernel, but part of the MPS build is host-side fp64 work]
    TN_CONTENTION=T(1.0, 1.2, 1.5),        # [assumption: 4 workers share one A100 + 16 cores; measured: 2 chi-64
                                           #  workers on one TITAN Xp showed no slowdown (tn_v2/timing_pair.jsonl)]
    A_STATE=T(2.0, 2.6, 3.2),              # [measured fits of the MPS-build time vs norb: 2.0 (chi 32, norb 15-18),
                                           #  2.2 (chi 64), 2.6 (chi 128), 3.2 (chi 256) over norb 15-29; central =
                                           #  the chi-128 fit (every size truncated, discarded weight 1e-4 to 4e-3);
                                           #  the chi-256 fit is inflated because the 15-17-orbital states are nearly
                                           #  exact at chi 256 (discarded weight ~1e-5); theory: n^2 at fixed chi]
    A_EXPECT=T(3.5, 4.0, 4.2),             # [measured fit norb 15-29: 3.6-4.2; QC-MPO bond ~n^2 over n sites]
    B_CHI_STATE=T(1.3, 1.8, 3.0),          # [measured chi exponents at n29: 1.2 (64->128), 1.8 (128->256), 1.3
                                           #  (256->512, one molecule at lower host load, probably low); used above
                                           #  chi 512 (sampling MPS only); 3.0 = the asymptotic chi^3]
    # ---- RL fine-tuning (GRPO, variational-energy reward) ----
    MOLS_PER_STEP=T(12, 16, 20),           # [measured 12 (grpo_n29_tn); assumption 16 at full scale]
    RL_STEPS_B1=T(90, 150, 225),           # [assumption] curriculum stage norb 15-19 (exact to 18, TN at 19)
    RL_STEPS_B2=T(180, 300, 450),          # [assumption] norb 20-29 (TN rewards); measured runs: 30-50 steps
    RL_STEPS_B3=T(270, 450, 675),          # [assumption] norb 30-37, 70% of the dataset
    RL_STEPS_B4=T(180, 300, 450),          # [assumption] norb 38-44
    RL_CHI_B2=T(64, 128, 256),             # [measured accuracy: dm engine median |err| 1.8 / 0.42 mHa at chi 64 / 128
    RL_CHI_B3=T(64, 128, 256),             #  (norb 15-18, tn_v2 acc/rank_dm.jsonl); at n29 chi 64 / 128 sit 2.1 / 0.54
    RL_CHI_B4=T(64, 128, 256),             #  mHa above chi 256 (n29_dm.jsonl, 11 states); chi 256 if larger molecules need it]
    RL_EVAL_FRAC=T(0.08, 0.15, 0.25),      # [measured: 7.6% (grpo_n29_tn), 21% (grpo_large_rl4) extra eval energies]
    RL_CAMPAIGNS=T(2, 3, 4),               # [assumption] main curriculum + seeds / reward-setting ablations
    RL_EFF=T(0.85, 0.75, 0.65),            # [assumption: synchronous GRPO steps wait for the slowest energy]
    EXACT_MULT=T(0.8, 1.0, 1.5),           # [measured: queue times 1.0-1.4x the benchmark; JIT on first call]
    HH_CAMPAIGNS=T(0.5, 1.0, 2.0),         # [assumption] heavy-hex variant for QPU circuits, norb <= 29 stage only
    N_DEV_RUNS=T(10, 20, 30),              # [assumption] short development runs at n29 scale (measured recipe)
    # ---- Evaluation: TN variational energies ----
    N_CAND=T(3, 4, 5),                     # [assumption] CCSD-truncated init, optimize=True label, NN, NN+RL (+1)
    EVAL_CHI=T(256, 256, 256),             # [measured: n29 chi 256 within ~0.1-0.4 mHa of convergence]
    CONV_FRAC=T(0.05, 0.10, 0.20),         # [assumption] fraction of test energies re-done at chi 512
    N_CKPT=T(20, 40, 80),                  # [assumption] checkpoints scored on a 100-molecule val subset at chi 128
    EVAL_EFF=T(0.9, 0.85, 0.8),            # [assumption] job packing / tails for queue-style evaluation jobs
    # ---- Evaluation: TN-sampled SQD (Lin et al. protocol) ----
    SQD_SUBSET_SCALE=T(0.6, 1.0, 1.5),     # [assumption] x the 300-molecule central subset
    SQD_CAND=T(3, 4, 5),                   # [assumption] parameter sets per molecule (5th: random / one-shot NN)
    SAMPLE_CHI_B2=T(256, 512, 1024),       # [measured: n29 discarded weight 9.8e-5 / 2.6e-5 at chi 256 / 512
    SAMPLE_CHI_B3=T(256, 512, 1024),       #  (tn_v2/n29_dm.jsonl); assumption: larger molecules need more]
    SAMPLE_CHI_B4=T(512, 1024, 2048),
    SHOTS=T(1e5, 1e5, 1e6),                # [paper: 1e5 samples used per state (1e6 cached for exact sampling)]
    SAMPLER_TFLOPS=T(10.0, 5.0, 2.0),      # [assumption: batched block-sparse GPU sampler (to be written), ~32 n
                                           #  chi^2 flops per shot]
    SQD_DIM=T(1e6, 4e6, 1.6e7),            # determinants per TN-sampled batch.  [paper: Lin et al.'s simulated
                                           #  subspaces stayed below the 4,000^2 = 1.6e7 cap: 242-929 strings per spin
                                           #  at the compressed init and <= 2,136 after TN optimization (1e5 exact
                                           #  samples; wan-hsuan-lucj/csv/improvement_quimb.xlsx).  assumption: our
                                           #  larger valence spaces (16-22 electrons per spin) give ~2,000 per spin;
                                           #  noisy or low-fidelity samples fill the cap (our MO-basis MPS samples at
                                           #  n29 were all unique: gen_logs/sqd_timing_n29.log)]
    SBD_DET_COST=T(8.0e-6, 1.24e-5, 1.93e-5),  # A100-s per determinant per diagonalization, GPU SBD solver
                                           #  [paper: Walkup et al. arXiv:2601.16169: N2 3.08e8 dets, 8x A100, 309 s
                                           #  -> 8.0e-6; H2O 6.28e8 dets, 32 Frontier nodes, 58.3 s -> 1.93e-5 per
                                           #  A100-equivalent; central = geometric mean]
    SBD_NLINK_EXP=T(1.0, 2.0, 2.5),        # growth of the per-determinant cost with nlink = nocc(nvirt+1) relative to
                                           #  the N2 anchor (110).  [assumption, central 2: reconciles our measured
                                           #  PySCF cost at norb 29 (nlink 224), the published GPU per-determinant
                                           #  costs at nlink ~100-110 and the published GPU/CPU speedup normalized to
                                           #  Perlmutter nodes (17-29x; implied here 27x); 1 = SBD's work grows more
                                           #  slowly than PySCF's link tables; 2.5 = steeper]
    SBD_EFF=T(0.8, 0.7, 0.6),              # [assumption: one node (4 A100) at d ~1e6-1e7 vs anchors at d ~3-6e8]
    # ---- QPU subset (QPU time is outside NERSC; the classical post-processing is GPU) ----
    QPU_MOLS=T(6, 9, 15),                  # [assumption] 2-5 test reactions x R/P/TS at norb <= 29 (<= 58 qubits)
    QPU_PARAM_SETS=T(2, 3, 4),             # [assumption] NN+RL (heavy-hex), compressed/optimize=True, random (+NN)
    QPU_RUNS=T(3, 5, 10),                  # [paper: 10 runs x 1e6 shots per circuit]
    QPU_SQD_RUNS=T(1, 3, 10),              # [paper: SQD per run (10); assumption: pooled (1) or 3 runs analysed]
    QPU_MIN_PER_1E6=T(4.5, 5.9, 6.5),      # [paper: IBM docs, 2 s + 0.00035 s/shot with default rep_delay]
    # ---- QSD ----
    QSD_FACTOR=T(1.0, 1.25, 1.5),          # [assumption: Fang Sun's estimate, QSD = 1.0-1.5x the SQD budget]
    # ---- Inference, timing, per-instance baselines, porting ----
    N_INF_CKPT=T(10, 30, 60),              # [assumption] checkpoints run over all 30,205 molecules
    INF_A100_S=T(0.01, 0.03, 0.1),         # [measured: 0.06-0.48 s per molecule on 4 CPU threads; GPU batched]
    TIMING_MOLS=T(100, 300, 600),          # [assumption] per-instance optimize=True timing on GPU
    TIMING_A100_S=T(30.0, 60.0, 120.0),    # [measured: 35-107 s per molecule on one CPU core at n29]
    BASE_MOLS=T(3, 6, 12),                 # [assumption] molecules for the per-instance TN-optimization baseline
    CALIB_NODE_H=T(20.0, 40.0, 80.0),      # [assumption] GPU SBD build/benchmark, dm-engine A100 scaling, training
    # ---- Contingency ----
    CONTINGENCY=T(0.12, 0.175, 0.22),      # [assumption] restarts, failed jobs, ablations; the n29 RL run needed
                                           #  ~1 node-h of failed attempts on top of 2.55 node-h
    # ---- CPU fallback (not in the GPU total) ----
    CPU_CORE_S_PER_DET_N29=T(0.029, 0.045, 0.063),  # core-s per determinant per diagonalization, PySCF selected CI
                                           #  [measured: gen_logs/sqd_timing_n29.log, 1 batch x 1 iteration, 4 threads,
                                           #  d = 0.9-1.5e5: 1,095-1,919 s]; scaled by nlink^2 (link tables) and
                                           #  linearly in d (extrapolated x100 to d = 1e7)
    # ---- Orbital-ordering labels (outline part 1) ----
    ORD_MOLS=T(0, 0, 0),           # not in this budget (Fang, 4 Oct 2026)
    _ORD_MOLS_V1=T(500, 1000, 2000),           # [assumption] molecules with ordering labels, drawn like the train pool
    ORD_LAYOUTS=T(10, 20, 50),             # [assumption] candidate layouts per molecule and topology (pilot: 15)
)
CPU_EFF = 0.7                # [assumption] parallel efficiency of PySCF/Dice on a 128-core CPU node
PROTO_STATES = 1000          # [assumption] states for CPU prototyping of SQD at reduced settings (CPU note only)
REQUEST_NODE_H = 10000      # request = central plan rounded to 10,000 (Fang, 4 Oct 2026)
REQUEST_PCTL = 0.65          # [assumption] request at about the 65th percentile of the Monte-Carlo total, because the
                             #  two largest unit costs (GPU SQD solver, TN engine on A100) are not yet benchmarked


def scenario(i):
    return {k: v[i] for k, v in PARAMS.items()}


# =================================================================================================================
# 3. Unit-cost functions
# =================================================================================================================
def _loglog(table, chi, hi_exp, lo_exp=1.0):
    """Piecewise power law through measured (chi, t) points; power-law extrapolation outside."""
    xs = sorted(table)
    if chi <= xs[0]:
        return table[xs[0]] * (chi / xs[0]) ** lo_exp
    if chi >= xs[-1]:
        return table[xs[-1]] * (chi / xs[-1]) ** hi_exp
    for a, b in zip(xs, xs[1:]):
        if a <= chi <= b:
            e = math.log(table[b] / table[a]) / math.log(b / a)
            return table[a] * (chi / a) ** e
    raise ValueError(chi)


def t_state_a100(n, chi, p):
    """Seconds of MPS build per state on an A100 for one worker (no co-scheduling)."""
    return _loglog(T_STATE_N29_TXP, chi, p["B_CHI_STATE"]) / p["A100_STATE_SPEEDUP"] * (n / 29) ** p["A_STATE"]


def t_expect_4thr(n, chi, p):
    return _loglog(T_EXPECT_N29_1THR, chi, 1.6) / B2_4THR_SPEEDUP * (n / 29) ** p["A_EXPECT"]


def energy_a100_s(n, chi, p, energies_per_mpo, exact_max=19):
    """A100-seconds per LUCJ variational energy <psi|H|psi>: exact on the GPU for norb <= exact_max, else TN dm
    engine with WORKERS_PER_GPU concurrent workers per A100 (4 block2 threads each)."""
    if n <= exact_max:
        return EXACT_A100_S[n] * p["EXACT_MULT"]
    worker = (t_state_a100(n, chi, p) + t_expect_4thr(n, chi, p)
              + T_MPO_N29 * (n / 29) ** A_MPO / energies_per_mpo)
    return worker * p["TN_CONTENTION"] / WORKERS_PER_GPU


def nlink(n):
    no = NOCC_MED[n]
    return no * (n - no + 1)


def diag_a100_s(n, d, p):
    """A100-seconds for one SQD subspace diagonalization of d determinants with the GPU SBD solver."""
    return p["SBD_DET_COST"] * (nlink(n) / NLINK_N2) ** p["SBD_NLINK_EXP"] * d / p["SBD_EFF"]


def diag_cpu_core_s(n, d, p):
    """CPU fallback: PySCF selected-CI core-seconds for one diagonalization."""
    return p["CPU_CORE_S_PER_DET_N29"] * (nlink(n) / nlink(29)) ** 2 * d


def gpu_vs_cpu_node_ratio(n, p):
    """Perlmutter CPU node-seconds / GPU node-seconds for the same diagonalization (the model's implied speedup)."""
    cpu = diag_cpu_core_s(n, 1.0, p) / (CORES_PER_CPU_NODE * CPU_EFF)
    gpu = diag_a100_s(n, 1.0, p) / A100_PER_NODE
    return cpu / gpu


def published_node_ratio():
    """Walkup et al.'s per-Frontier-node speedup normalized to a Perlmutter GPU node vs a Perlmutter CPU node."""
    gpu_frac = A100_PER_NODE * A100_IN_GCD / 8.0                    # Perlmutter GPU node / Frontier GPU side
    lo = WALKUP_SPEEDUP_PER_FRONTIER_NODE[0] * gpu_frac / (CORES_PER_CPU_NODE / FRONTIER_CPU_CORES[0])
    hi = WALKUP_SPEEDUP_PER_FRONTIER_NODE[1] * gpu_frac / (CORES_PER_CPU_NODE / FRONTIER_CPU_CORES[1])
    return lo, hi


def sample_chi(n, p):
    if n <= 29:
        return p["SAMPLE_CHI_B2"]
    return p["SAMPLE_CHI_B3"] if n <= 37 else p["SAMPLE_CHI_B4"]


def sqd_state_parts(n, p, d=None, n_diag=N_BATCH * TN_SQD_ITERS, shots=None, chi=None):
    """A100-seconds for one SQD-evaluated state: (state preparation, sampling, diagonalizations)."""
    d = p["SQD_DIM"] if d is None else d
    shots = p["SHOTS"] if shots is None else shots
    chi = sample_chi(n, p) if chi is None else chi
    if n <= 19:                                     # exact state vector on the GPU
        prep, samp = EXACT_A100_S[n] * p["EXACT_MULT"], 0.0
    else:                                           # one MPS per GPU at sampling chi (GPU-bound, no co-scheduling)
        prep = t_state_a100(n, chi, p)
        samp = shots * 32.0 * n * chi ** 2 / (p["SAMPLER_TFLOPS"] * 1e12)
    return prep, samp, n_diag * diag_a100_s(n, d, p)


def base_eval_a100_s(p, n_diag=BASE_BATCHES, d=BASE_D):
    """One objective evaluation of the per-molecule TN optimization (MPS chi 128 + 1e4 shots + n_diag diags)."""
    return (t_state_a100(BASE_NORB, BASE_CHI, p) * p["TN_CONTENTION"] / WORKERS_PER_GPU
            + BASE_SHOTS * 32.0 * BASE_NORB * BASE_CHI ** 2 / (p["SAMPLER_TFLOPS"] * 1e12)
            + n_diag * diag_a100_s(BASE_NORB, d, p))


def wmean(dist, f):
    tot = sum(dist.values())
    return sum(c * f(n) for n, c in dist.items()) / tot


def bucket_dist(dist, b):
    lo, hi = BUCKETS[b]
    return {n: c for n, c in dist.items() if lo <= n <= hi and c > 0}


# =================================================================================================================
# 4. Budget (A100-hours per line item)
# =================================================================================================================
def budget(p):
    L = {}        # line item -> A100-hours
    X = {}        # extra outputs (counts, QPU minutes, CPU fallback)
    h = 1 / 3600.0

    # ---- Pretraining ----
    epoch_a100_h = p["EPOCH_LAB_GPU_H"] / p["LAB_TO_A100"]
    L["pre_prod"] = p["N_PROD"] * p["EPOCHS_PROD"] * p["MODEL_MULT_PROD"] * epoch_a100_h
    L["pre_sweep"] = p["N_SWEEP"] * p["EPOCHS_SWEEP"] * p["MODEL_MULT_SWEEP"] * epoch_a100_h
    X["epoch_a100_h"] = epoch_a100_h

    # ---- Orbital-ordering labels (MPS/LUCJ energy judge; exact <= 18 orbitals, TN chi 128 above) ----
    k_ord = p["ORD_LAYOUTS"] * ORD_TOPOLOGIES
    ord_unit = wmean(TRAIN, lambda n: energy_a100_s(n, 128, p, k_ord, EXACT_MAX_NORB_RL))
    X["ord_energies"] = p["ORD_MOLS"] * k_ord
    X["ord_unit"] = ord_unit
    L["ord_labels"] = X["ord_energies"] * ord_unit * h / p["EVAL_EFF"]
    # CPU part: one optimize=True compressed-DF fit per layout, 35/70/107 core-s at n29, x (n/29)^4
    X["ord_cpu_core_h"] = X["ord_energies"] * wmean(TRAIN, lambda n: (n / 29) ** 4) * 70.0 / 3600

    # ---- RL fine-tuning ----
    e_per_step = p["MOLS_PER_STEP"] * GROUP * (1 + p["RL_EVAL_FRAC"])
    chi_of = {"B1": 128, "B2": p["RL_CHI_B2"], "B3": p["RL_CHI_B3"], "B4": p["RL_CHI_B4"]}
    rl_energies = 0.0
    for b in ("B1", "B2", "B3", "B4"):
        pool = bucket_dist(TRAIN, b)
        unit = wmean(pool, lambda n: energy_a100_s(n, chi_of[b], p, GROUP, EXACT_MAX_NORB_RL))
        E = p[f"RL_STEPS_{b}"] * e_per_step
        X[f"rl_unit_{b}"] = unit
        X[f"rl_energies_{b}"] = E * p["RL_CAMPAIGNS"]
        L[f"rl_{b}"] = E * unit * h / p["RL_EFF"] * p["RL_CAMPAIGNS"]
        rl_energies += E * p["RL_CAMPAIGNS"]
        if b == "B2":
            L["rl_heavyhex"] = p["HH_CAMPAIGNS"] * E * unit * h / p["RL_EFF"]
            X["rl_energies_hh"] = p["HH_CAMPAIGNS"] * E
            rl_energies += p["HH_CAMPAIGNS"] * E
    dev_E = p["N_DEV_RUNS"] * DEV_RUN_ENERGIES
    X["rl_unit_dev"] = energy_a100_s(29, 128, p, GROUP)
    L["rl_dev"] = dev_E * X["rl_unit_dev"] * h / p["RL_EFF"]
    X["rl_energies"] = rl_energies + dev_E
    # CPU fallback for the TN rewards: v1 zip-up chi 64 on CPU nodes, scaled with the measured chi-64 total exponent
    cpu_rl = 0.0
    for b in ("B2", "B3", "B4"):
        pool = bucket_dist(TRAIN, b)
        per_e = wmean(pool, lambda n: (n / 29) ** 3.1) / V1_CHI64_CPU_ENERGIES_PER_NODE_H
        cpu_rl += (X[f"rl_energies_{b}"] + (X["rl_energies_hh"] if b == "B2" else 0.0)) * per_e
    X["rl_cpu_fallback_node_h"] = cpu_rl + dev_E / V1_CHI64_CPU_ENERGIES_PER_NODE_H

    # ---- Evaluation: TN variational energies ----
    ncand = p["N_CAND"]
    test_e = sum(c * ncand * energy_a100_s(n, p["EVAL_CHI"], p, ncand) for n, c in TEST.items())
    val_e = sum(c * ncand * energy_a100_s(n, p["EVAL_CHI"], p, ncand) for n, c in VAL.items())
    nval = sum(VAL.values())
    sel_e = p["N_CKPT"] * sum((100.0 * c / nval) * energy_a100_s(n, 128, p, 1) for n, c in VAL.items())
    conv_e = p["CONV_FRAC"] * sum(c * ncand * energy_a100_s(n, 512, p, ncand) for n, c in TEST.items())
    L["eval_tn_energy"] = (test_e + val_e + sel_e + conv_e) * h / p["EVAL_EFF"]
    X["eval_energies"] = ((sum(TEST.values()) + nval) * ncand + p["N_CKPT"] * 100
                          + p["CONV_FRAC"] * sum(TEST.values()) * ncand)

    # ---- Evaluation: TN-sampled SQD (Lin et al. protocol, adapted) ----
    prep = samp = diag = 0.0
    cpu_core_s = cpu_core_s_red = 0.0
    n_states = 0.0
    for n, c in SQD_SUBSET.items():
        k = c * p["SQD_SUBSET_SCALE"] * p["SQD_CAND"]
        a, s, dg = sqd_state_parts(n, p)
        prep += k * a
        samp += k * s
        diag += k * dg
        cpu_core_s += k * N_BATCH * TN_SQD_ITERS * diag_cpu_core_s(n, p["SQD_DIM"], p)
        cpu_core_s_red += k * N_BATCH * TN_SQD_ITERS * diag_cpu_core_s(n, min(p["SQD_DIM"], 4e6), p)
        n_states += k
    L["sqd_tn_prep"] = prep * h
    L["sqd_tn_sample"] = samp * h
    L["sqd_tn_diag"] = diag * h
    X["sqd_states"] = n_states
    X["n_diag_tn"] = n_states * N_BATCH * TN_SQD_ITERS

    # ---- Evaluation: QPU subset, classical part (GPU) ----
    circuits = p["QPU_MOLS"] * p["QPU_PARAM_SETS"]
    sqd_runs = min(p["QPU_SQD_RUNS"], p["QPU_RUNS"])
    n_diag_circ = sqd_runs * QPU_CR_ITERS * N_BATCH
    L["sqd_qpu_post"] = circuits * n_diag_circ * diag_a100_s(QPU_NORB, QPU_D, p) * h
    a, s, dg = sqd_state_parts(QPU_NORB, p, chi=p["SAMPLE_CHI_B2"])        # noiseless TN reference, same circuits
    L["sqd_qpu_emul"] = circuits * (a + s + dg) * h
    for dq, red in ((QPU_D, False), (min(QPU_D, 4e6), True)):
        v = (circuits * n_diag_circ * diag_cpu_core_s(QPU_NORB, dq, p)
             + circuits * N_BATCH * diag_cpu_core_s(QPU_NORB, min(p["SQD_DIM"], dq), p))
        if red:
            cpu_core_s_red += v
        else:
            cpu_core_s += v
    X["qpu_circuits"] = circuits
    X["n_diag_qpu"] = circuits * (n_diag_circ + N_BATCH)
    X["qpu_shots"] = circuits * p["QPU_RUNS"] * 1e6
    X["qpu_min_sqd"] = circuits * p["QPU_RUNS"] * p["QPU_MIN_PER_1E6"]
    X["qpu_min_qsd"] = X["qpu_min_sqd"] * p["QSD_FACTOR"]

    sqd_total = sum(L[k] for k in ("sqd_tn_prep", "sqd_tn_sample", "sqd_tn_diag", "sqd_qpu_post", "sqd_qpu_emul"))
    L["qsd"] = p["QSD_FACTOR"] * sqd_total
    X["sqd_cpu_fallback_node_h"] = cpu_core_s / 3600 / CORES_PER_CPU_NODE / CPU_EFF
    X["qsd_cpu_fallback_node_h"] = X["sqd_cpu_fallback_node_h"] * p["QSD_FACTOR"]
    X["sqd_cpu_fallback_red_node_h"] = cpu_core_s_red / 3600 / CORES_PER_CPU_NODE / CPU_EFF
    X["sqd_gpu_diag_node_h"] = (L["sqd_tn_diag"] + L["sqd_qpu_post"]
                                + circuits * dg * h) / A100_PER_NODE

    # ---- Inference, timing, per-instance baselines ----
    inf = p["N_INF_CKPT"] * N_ALL * p["INF_A100_S"] * h
    timing = p["TIMING_MOLS"] * p["TIMING_A100_S"] * h
    base = p["BASE_MOLS"] * BASE_EVALS * base_eval_a100_s(p) * h
    L["inf_timing_base"] = inf + timing + base
    X["base_node_h_per_mol"] = BASE_EVALS * base_eval_a100_s(p) * h / A100_PER_NODE
    # range for the per-molecule TN-optimization cost (Lin et al.'s NOMAD loop, 500 evaluations):
    #   Lin-like: MPS-sampled subspaces of 101-543 strings per spin with identical batches -> 1 diag at d = 3e5
    #   ours at the central subspace: 10 distinct batches at d = 4e6;  full cap: 10 batches at d = 1.6e7
    X["base_lin_node_h"] = BASE_EVALS * base_eval_a100_s(p, 1, 3e5) * h / A100_PER_NODE
    X["base_mid_node_h"] = BASE_EVALS * base_eval_a100_s(p, N_BATCH, 4e6) * h / A100_PER_NODE
    X["base_full_node_h_per_mol"] = BASE_EVALS * base_eval_a100_s(p, N_BATCH, 1.6e7) * h / A100_PER_NODE
    L["calib"] = p["CALIB_NODE_H"] * A100_PER_NODE

    sub = sum(L.values())
    L["contingency"] = p["CONTINGENCY"] * sub
    return L, X


GROUPS = [  # (label, keys) for the summary table
    ("Pretraining: production runs", ["pre_prod"]),
    ("Pretraining: architecture/hyperparameter sweeps", ["pre_sweep"]),
    ("RL: norb 15-19 (exact to 18, TN at 19)", ["rl_B1"]),
    ("RL: TN rewards, norb 20-29", ["rl_B2"]),
    ("RL: TN rewards, norb 30-37", ["rl_B3"]),
    ("RL: TN rewards, norb 38-44", ["rl_B4"]),
    ("RL: heavy-hex variant (QPU circuits)", ["rl_heavyhex"]),
    ("RL: development runs", ["rl_dev"]),
    ("Eval: TN variational energies", ["eval_tn_energy"]),
    ("Eval: SQD, TN-sampled (MPS + samples)", ["sqd_tn_prep", "sqd_tn_sample"]),
    ("Eval: SQD, TN-sampled (diagonalization)", ["sqd_tn_diag"]),
    ("Eval: SQD on QPU samples + TN emulation", ["sqd_qpu_post", "sqd_qpu_emul"]),
    ("Eval: QSD (= factor x SQD)", ["qsd"]),
    ("Inference, timing, per-instance baselines", ["inf_timing_base"]),
    ("Porting & calibration on Perlmutter", ["calib"]),
    ("Contingency (restarts, failures, ablations)", ["contingency"]),
]
SECTIONS = [  # (label, keys) for the 'Planned GPU usage' grouping
    ("Pretraining + sweeps", ["pre_prod", "pre_sweep"]),
    ("RL fine-tuning (TN/exact rewards)", ["rl_B1", "rl_B2", "rl_B3", "rl_B4", "rl_heavyhex", "rl_dev"]),
    ("TN energies (test/val/model selection)", ["eval_tn_energy"]),
    ("SQD testing (TN-sampled + QPU post-processing)", ["sqd_tn_prep", "sqd_tn_sample", "sqd_tn_diag",
                                                       "sqd_qpu_post", "sqd_qpu_emul"]),
    ("QSD testing", ["qsd"]),
    ("Inference/timing/baselines + porting", ["inf_timing_base", "calib"]),
    ("Contingency", ["contingency"]),
]
OLD_REQUEST = [  # the earlier proposal text (10,000 node-h) mapped onto the v2 line items
    ("Supervised pretraining on teacher LUCJ parameters", 2500, ["pre_prod"]),
    ("Architecture and hyperparameter sweeps", 2000, ["pre_sweep"]),
    ("RL fine-tuning with QSCI energy as reward", 3000, ["rl_B1", "rl_B2", "rl_B3", "rl_B4", "rl_heavyhex", "rl_dev"]),
    ("Inference sweeps, generalization benchmarks, timing", 1000, ["inf_timing_base", "calib", "eval_tn_energy"]),
    ("(not in the old text) SQD and QSD tests", 0, ["sqd_tn_prep", "sqd_tn_sample", "sqd_tn_diag", "sqd_qpu_post",
                                                    "sqd_qpu_emul", "qsd"]),
    ("Restarts, failure recovery, ablations, production", 1500, ["contingency"]),
]


def node_h(L, keys):
    return sum(L[k] for k in keys) / A100_PER_NODE


# =================================================================================================================
# 5. Monte Carlo and sensitivity
# =================================================================================================================
def sample_params(rng):
    p = {}
    for k, (lo, c, hi) in PARAMS.items():
        a, b = min(lo, hi), max(lo, hi)
        if a == b:
            p[k] = c
        elif a <= 0:
            p[k] = rng.triangular(a, b, c)
        else:
            p[k] = math.exp(rng.triangular(math.log(a), math.log(b), math.log(c)))
    return p


def total_node_h(p):
    L, _ = budget(p)
    return sum(L.values()) / A100_PER_NODE


def formulas(p, L, X):
    """Central-scenario formula strings with explicit unit conversions (/3600 s per h, /4 A100 per node)."""
    e_step = p["MOLS_PER_STEP"] * GROUP * (1 + p["RL_EVAL_FRAC"])
    f = {}
    f["pre_prod"] = (f"{p['N_PROD']:.0f} runs × {p['EPOCHS_PROD']:.0f} epochs × {p['MODEL_MULT_PROD']} (model size) "
                     f"× {X['epoch_a100_h']:.2f} A100-h ÷ 4")
    f["pre_sweep"] = (f"{p['N_SWEEP']:.0f} configs × {p['EPOCHS_SWEEP']:.0f} epochs × {p['MODEL_MULT_SWEEP']} "
                      f"× {X['epoch_a100_h']:.2f} A100-h ÷ 4")
    f["ord_labels"] = (f"{p['ORD_MOLS']:,.0f} molecules × {p['ORD_LAYOUTS']:.0f} layouts × {ORD_TOPOLOGIES} topologies "
                       f"= {X['ord_energies']:,.0f} energies × {X['ord_unit']:.1f} A100-s ÷ {p['EVAL_EFF']} "
                       f"÷ 3600 ÷ 4")
    for b in ("B1", "B2", "B3", "B4"):
        f[f"rl_{b}"] = (f"{p['RL_CAMPAIGNS']:.0f} × {p[f'RL_STEPS_{b}']:.0f} steps × {p['MOLS_PER_STEP']:.0f} mol × "
                        f"{GROUP} × {1 + p['RL_EVAL_FRAC']:.2f} = {X[f'rl_energies_{b}']:,.0f} energies × "
                        f"{X[f'rl_unit_{b}']:.1f} A100-s ÷ {p['RL_EFF']} ÷ 3600 ÷ 4")
    f["rl_heavyhex"] = (f"{p['HH_CAMPAIGNS']:.0f} × the 20–29-orbital stage ({p['RL_STEPS_B2'] * e_step:,.0f} energies "
                        f"× {X['rl_unit_B2']:.1f} A100-s ÷ {p['RL_EFF']} ÷ 3600 ÷ 4)")
    f["rl_dev"] = (f"{p['N_DEV_RUNS']:.0f} runs × {DEV_RUN_ENERGIES:,} energies × {X['rl_unit_dev']:.1f} A100-s "
                   f"÷ {p['RL_EFF']} ÷ 3600 ÷ 4")
    f["eval_tn_energy"] = (f"(1,174 test + 675 val) × {p['N_CAND']:.0f} at χ{p['EVAL_CHI']:.0f} + {p['N_CKPT']:.0f} "
                           f"checkpoints × 100 val at χ128 + {p['CONV_FRAC']:.0%} of test energies at χ512; "
                           f"÷ {p['EVAL_EFF']}")
    f["sqd_tn"] = (f"{X['sqd_states']:,.0f} states ({sum(SQD_SUBSET.values()) * p['SQD_SUBSET_SCALE']:.0f} test "
                   f"molecules × {p['SQD_CAND']:.0f} parameter sets) × [MPS at χ{p['SAMPLE_CHI_B2']:.0f}–"
                   f"{p['SAMPLE_CHI_B4']:.0f} + {p['SHOTS']:.0e} shots + {N_BATCH} diagonalizations at d = "
                   f"{p['SQD_DIM']:.1e}]")
    f["sqd_qpu"] = (f"{X['qpu_circuits']:.0f} circuits × {min(p['QPU_SQD_RUNS'], p['QPU_RUNS']):.0f} SQD runs × "
                    f"{QPU_CR_ITERS} recovery iterations × {N_BATCH} batches at d = {QPU_D:.1e} (norb {QPU_NORB}) "
                    f"+ noiseless TN emulation of the same circuits")
    f["qsd"] = f"{p['QSD_FACTOR']} × the three SQD lines (Fang's estimate 1.0–1.5×)"
    f["inf_timing_base"] = (f"inference {p['N_INF_CKPT']:.0f} ckpt × {N_ALL:,} mol × {p['INF_A100_S']} A100-s; "
                            f"timing {p['TIMING_MOLS']:.0f} × {p['TIMING_A100_S']:.0f} A100-s; NOMAD baseline "
                            f"{p['BASE_MOLS']:.0f} mol × {X['base_node_h_per_mol']:.1f} node-h")
    f["calib"] = f"{p['CALIB_NODE_H']:.0f} node-h"
    f["contingency"] = f"{p['CONTINGENCY']:.0%} of the subtotal"
    return f


def pct(vals, q):
    s = sorted(vals)
    return s[int(q * (len(s) - 1))]


def main(md=False, n_mc=4000):
    S = [scenario(i) for i in range(3)]
    R = [budget(p) for p in S]
    Lc, Xc = R[1]

    print("=" * 112)
    print("ML-LUCJ GPU budget model v2.  Unit: Perlmutter GPU node-hour (4x A100).  Central = all inputs at central "
          "values;")
    print("Low/High = Monte-Carlo P10/P90 (inputs drawn log-triangular between their low and high ends, mode = central);")
    print("corners = every input at its cheap / expensive end simultaneously (shown for transparency only).")
    print("=" * 112)
    print(f"dataset: {N_ALL} molecules; train {sum(TRAIN.values())}, val {sum(VAL.values())}, "
          f"test {sum(TEST.values())}; central SQD subset {sum(SQD_SUBSET.values())} molecules")

    # ---- unit costs ----
    print("\nUnit cost: A100-seconds per LUCJ variational energy (exact norb<=19; TN dm engine, 4 workers/A100); "
          "corner low/central/high")
    print(f"{'norb':>5} {'qubits':>6} | {'chi128':>16} | {'chi256':>16} | {'chi512 cen':>10}")
    for n in (17, 18, 19, 20, 25, 29, 33, 37, 40, 44):
        r1 = "/".join(f"{energy_a100_s(n, 128, p, GROUP):.0f}" for p in S)
        r2 = "/".join(f"{energy_a100_s(n, 256, p, 4):.0f}" for p in S)
        r3 = f"{energy_a100_s(n, 512, S[1], 4):.0f}"
        print(f"{n:>5} {2 * n:>6} | {r1:>16} | {r2:>16} | {r3:>10}")
    print("  check vs measured: v1 zip-up chi256 at n29 on A100, 5 workers/GPU: 59-80 A100-s per energy "
          "(pretrain/opt_true/results/energy_n29*_chi256.json); model chi256 n29 central "
          f"{energy_a100_s(29, 256, S[1], 4):.0f}")
    c = S[1]
    ts, te, tm = t_state_a100(29, 128, c), t_expect_4thr(29, 128, c), T_MPO_N29 / GROUP
    print(f"  TN worker wall time at n29 chi128, central: state {ts:.2f} s + <H> (4 thr) {te:.2f} s + MPO/8 {tm:.2f} s"
          f" = {ts + te + tm:.1f} s x {c['TN_CONTENTION']} / {WORKERS_PER_GPU} = "
          f"{(ts + te + tm) * c['TN_CONTENTION'] / WORKERS_PER_GPU:.1f} A100-s")
    print(f"  formula: cost(n) = [{ts:.1f}*(n/29)^{c['A_STATE']} + {te:.1f}*(n/29)^{c['A_EXPECT']} + "
          f"{T_MPO_N29:.0f}*(n/29)^{A_MPO:.0f}/8] x {c['TN_CONTENTION']}/{WORKERS_PER_GPU}  ->  " + ", ".join(
              f"n{n}: {energy_a100_s(n, 128, c, GROUP):.1f}" for n in (20, 25, 29, 33, 37, 40, 44)))
    print("\nUnit cost: one SQD-evaluated state (1e5 shots, 10 batches, QSCI), A100-hours (corner low/central/high)")
    for n in (17, 25, 29, 33, 37, 40, 44):
        parts = [sqd_state_parts(n, p) for p in S]
        tot = "/".join(f"{sum(x) / 3600:.2f}" for x in parts)
        a, s, dg = parts[1]
        print(f"  norb {n} (chi {sample_chi(n, S[1])}, nlink {nlink(n)}): {tot}; central = MPS {a / 3600:.2f} + "
              f"samples {s / 3600:.3f} + diag {dg / 3600:.2f};  1 diag at d=1e7: {diag_a100_s(n, 1e7, S[1]) / 60:.1f}"
              f" A100-min (corners {diag_a100_s(n, 1e7, S[0]) / 60:.1f}-{diag_a100_s(n, 1e7, S[2]) / 60:.1f}); "
              f"at d=4e6: {diag_a100_s(n, 4e6, S[1]) / 60:.1f}")
    print("  CPU fallback (PySCF), one diagonalization at d=1e7: " + ", ".join(
        f"norb {n}: {diag_cpu_core_s(n, 1e7, S[1]) / 3600:.0f} core-h" for n in (29, 33, 39, 44)))
    plo, phi = published_node_ratio()
    print(f"  consistency: implied Perlmutter CPU-node-s / GPU-node-s for the same diagonalization (central): " +
          ", ".join(f"norb {n}: {gpu_vs_cpu_node_ratio(n, S[1]):.0f}x" for n in (17, 29, 33, 44)) +
          f";  published (Walkup, normalized to Perlmutter nodes): {plo:.0f}-{phi:.0f}x")
    for e in (1.0, 2.0):
        q = dict(S[1]); q["SBD_NLINK_EXP"] = e
        print(f"    with SBD_NLINK_EXP={e}: " + ", ".join(f"norb {n}: {gpu_vs_cpu_node_ratio(n, q):.0f}x"
                                                         for n in (17, 29, 33, 44)))
    print(f"  SQD state at n29 (central) / one TN reward energy at n29 chi128: "
          f"{sum(sqd_state_parts(29, S[1])) / energy_a100_s(29, 128, S[1], GROUP):.0f}x")

    # ---- RL detail ----
    print("\nRL detail: mean A100-s per reward energy per curriculum stage (corner low/central/high); central energies")
    for b in ("B1", "B2", "B3", "B4"):
        print(f"  {b} norb {BUCKETS[b][0]}-{BUCKETS[b][1]} (train pool {sum(bucket_dist(TRAIN, b).values()):,}): "
              f"{R[0][1][f'rl_unit_{b}']:.1f} / {Xc[f'rl_unit_{b}']:.1f} / {R[2][1][f'rl_unit_{b}']:.1f} A100-s; "
              f"energies {Xc[f'rl_energies_{b}']:,.0f}")
    print(f"  total RL energies (central, incl. heavy-hex + development runs): {Xc['rl_energies']:,.0f}")
    print(f"  ordering labels: {Xc['ord_energies']:,.0f} energies x {Xc['ord_unit']:.1f} A100-s (train-pool mix)")

    # ---- Monte Carlo ----
    rng = random.Random(20261004)
    draws = []
    for _ in range(n_mc):
        L, X = budget(sample_params(rng))
        draws.append((L, X))
    tot_mc = [sum(L.values()) / A100_PER_NODE for L, _ in draws]

    def line_mc(keys):
        v = [sum(L[k] for k in keys) / A100_PER_NODE for L, _ in draws]
        return pct(v, 0.10), pct(v, 0.50), pct(v, 0.90)

    # ---- formulas ----
    F = formulas(S[1], Lc, Xc)
    print("\nCentral formulas:")
    for k, s in F.items():
        print(f"  {k:16s} {s}")

    # ---- line items ----
    print("\nLine items, node-hours: central | MC P10 / P50 / P90 | corners low / high | central A100-h")
    print(f"{'item':48s} | {'central':>7} | {'P10':>6} {'P50':>6} {'P90':>6} | {'c.low':>6} {'c.high':>7} | "
          f"{'A100-h':>7}")
    for label, keys in GROUPS:
        c_ = node_h(Lc, keys)
        lo, med, hi = line_mc(keys)
        cl, ch = node_h(R[0][0], keys), node_h(R[2][0], keys)
        print(f"{label:48s} | {c_:7.0f} | {lo:6.0f} {med:6.0f} {hi:6.0f} | {cl:6.0f} {ch:7.0f} | {c_ * 4:7.0f}")
    tot = [sum(R[i][0].values()) / A100_PER_NODE for i in range(3)]
    print(f"{'TOTAL':48s} | {tot[1]:7.0f} | {pct(tot_mc, .1):6.0f} {pct(tot_mc, .5):6.0f} {pct(tot_mc, .9):6.0f} | "
          f"{tot[0]:6.0f} {tot[2]:7.0f} | {tot[1] * 4:7.0f}")
    print(f"  MC mean {statistics.mean(tot_mc):.0f}; P05 {pct(tot_mc, .05):.0f}; P95 {pct(tot_mc, .95):.0f}; "
          f"P99 {pct(tot_mc, .99):.0f}; P(total > 10,000) = {sum(t > 10000 for t in tot_mc) / len(tot_mc):.1%}; "
          f"P(total < central) = {sum(t < tot[1] for t in tot_mc) / len(tot_mc):.0%}")
    print("  MC percentiles: " + "; ".join(f"P{int(q * 100)} {pct(tot_mc, q):,.0f}"
                                           for q in (0.25, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8)))
    req = REQUEST_NODE_H
    print(f"  request: {req:,} node-h (central rounded); "
          f"P(total > request) = {sum(t > req for t in tot_mc) / len(tot_mc):.0%}; "
          f"P(total > 10,000) = {sum(t > 10000 for t in tot_mc) / len(tot_mc):.1%}")

    print("\nBy section (node-hours): central (share) | MC P10-P90")
    for label, keys in SECTIONS:
        c_ = node_h(Lc, keys)
        lo, _, hi = line_mc(keys)
        print(f"  {label:48s} {c_:7.0f} ({c_ / tot[1] * 100:4.1f}%) | {lo:6.0f}-{hi:6.0f}")
    print("\nOld request (10,000 node-h) vs this model (central node-h):")
    for label, old, keys in OLD_REQUEST:
        print(f"  {label:58s} old {old:6,}  new {node_h(Lc, keys):6,.0f}")

    # ---- sensitivity ----
    sens = []
    for k, (lo, c_, hi) in PARAMS.items():
        if lo == hi == c_ or k.startswith("CPU_"):
            continue
        lo_p, hi_p = dict(S[1]), dict(S[1])
        lo_p[k], hi_p[k] = lo, hi
        a, b = total_node_h(lo_p), total_node_h(hi_p)
        sens.append((b - a, k, a, b))
    sens.sort(reverse=True)
    print(f"\nOne-at-a-time sensitivity around central ({tot[1]:.0f} node-h): top drivers")
    for sw, k, a, b in sens[:14]:
        print(f"  {k:24s} {str(PARAMS[k]):34s} -> {a:6.0f} .. {b:6.0f} node-h (swing {sw:5.0f})")
    for k, v in (("SQD_DIM", 1.6e7), ("SBD_NLINK_EXP", 1.0), ("QSD_FACTOR", 1.5), ("QSD_FACTOR", 1.0),
                 ("B_CHI_STATE", 3.0)):
        q = dict(S[1]); q[k] = v
        print(f"  what-if {k}={v}: total {total_node_h(q):,.0f}")
    q = dict(S[1]); q["SQD_DIM"] = 1e7; q["SBD_NLINK_EXP"] = 1.0; q["A_STATE"] = 2.5
    print(f"  check: v1 central inputs, without the ordering line, reproduce the draft total 5,173: "
          f"{total_node_h(q) - node_h(budget(q)[0], ['ord_labels']) * (1 + q['CONTINGENCY']):,.0f}")
    for k, v in (("SQD_DIM", 4e6), ("SBD_NLINK_EXP", 2.0), ("A_STATE", 2.6)):
        q[k] = v
        print(f"    + {k}={v}: {total_node_h(q) - node_h(budget(q)[0], ['ord_labels']) * (1 + q['CONTINGENCY']):,.0f}")
    print(f"    + ordering-label line (with contingency): {total_node_h(q):,.0f}")

    # ---- outside the GPU request ----
    mcx = lambda key: (pct([x[key] for _, x in draws], .1), pct([x[key] for _, x in draws], .9))
    print("\nOutside the GPU request: central (MC P10-P90)")
    print(f"  QPU circuits {Xc['qpu_circuits']:.0f}; shots {Xc['qpu_shots']:.2e}")
    a, b = mcx("qpu_min_sqd"), mcx("qpu_min_qsd")
    qmin = Xc['qpu_min_sqd'] + Xc['qpu_min_qsd']
    print(f"  QPU minutes: SQD {Xc['qpu_min_sqd']:.0f} ({a[0]:.0f}-{a[1]:.0f}); QSD {Xc['qpu_min_qsd']:.0f} "
          f"({b[0]:.0f}-{b[1]:.0f}); total {qmin:.0f} min = {qmin / 60:.1f} h; list price ${qmin * 48 / 1e3:,.0f}k"
          f"-${qmin * 96 / 1e3:,.0f}k at $48-96/min; {qmin / 400:.1f} default OLCF QCUP windows (400 min / 28 days)")
    a, b = mcx("sqd_cpu_fallback_node_h"), mcx("qsd_cpu_fallback_node_h")
    print(f"  CPU fallback if the SQD diagonalizations run on CPU (PySCF): SQD {Xc['sqd_cpu_fallback_node_h']:,.0f} "
          f"({a[0]:,.0f}-{a[1]:,.0f}) CPU node-h; QSD {Xc['qsd_cpu_fallback_node_h']:,.0f} ({b[0]:,.0f}-{b[1]:,.0f})")
    a = mcx("sqd_gpu_diag_node_h")
    print(f"     ... this replaces {Xc['sqd_gpu_diag_node_h']:,.0f} ({a[0]:,.0f}-{a[1]:,.0f}) GPU node-h of SQD "
          f"diagonalization (x(1+QSD factor) incl. QSD = {Xc['sqd_gpu_diag_node_h'] * (1 + S[1]['QSD_FACTOR']):,.0f})")
    red = Xc["sqd_cpu_fallback_red_node_h"]
    print(f"     with max_dim capped at 2,000 (d <= 4e6, both tiers): SQD {red:,.0f} CPU node-h "
          f"({Xc['sqd_cpu_fallback_node_h'] / red:.1f}x cheaper); incl. QSD {red * (1 + S[1]['QSD_FACTOR']):,.0f}")
    a = mcx("rl_cpu_fallback_node_h")
    print(f"  CPU fallback for the RL TN rewards (v1 zip-up chi 64 on 128-core CPU nodes): "
          f"{Xc['rl_cpu_fallback_node_h']:,.0f} ({a[0]:,.0f}-{a[1]:,.0f}) CPU node-h")
    print(f"  per-instance TN-optimization baseline (Lin et al. NOMAD, 500 evaluations), node-h per molecule at "
          f"norb 29: budgeted (3 batches, d=1e6) {Xc['base_node_h_per_mol']:.1f}; Lin-like (1 distinct diag, d=3e5) "
          f"{Xc['base_lin_node_h']:.1f}; 10 diags at d=4e6 {Xc['base_mid_node_h']:.0f}; 10 at the cap 1.6e7 "
          f"{Xc['base_full_node_h_per_mol']:.0f}")
    print(f"  SQD diagonalizations (central): TN tier {Xc['n_diag_tn']:,.0f}; QPU tier {Xc['n_diag_qpu']:,.0f}; "
          f"total {Xc['n_diag_tn'] + Xc['n_diag_qpu']:,.0f}")
    # CPU prototyping at reduced settings: PROTO_STATES states x 3 batches at d = 1e6 (max_dim 1000), test-subset mix
    proto = [PROTO_STATES * wmean(SQD_SUBSET, lambda n: 3 * diag_cpu_core_s(n, 1e6, p)) / 3600 / CORES_PER_CPU_NODE
             / CPU_EFF for p in S]
    print(f"  CPU prototyping of SQD at reduced settings ({PROTO_STATES} states x 3 batches, d = 1e6): "
          + " / ".join(f"{v:,.0f}" for v in proto) + " CPU node-h (corner low/central/high)")
    # labels for val+test: optimize=True, 3 configs, 35-107 core-s at n29 scaled ~n^4
    lab = [(sum(TEST.values()) + sum(VAL.values())) * 3 * s * wmean({**TEST, **VAL}, lambda n: (n / 29) ** 4) / 3600
           for s in (35.0, 70.0, 107.0)]
    print("  optimize=True labels for val+test (3 configs, 35/70/107 core-s at n29, x(n/29)^4): "
          + " / ".join(f"{v:,.0f}" for v in lab) + " core-h")
    print(f"  optimize=True fits for the ordering labels (70 core-s at n29, x(n/29)^4): {Xc['ord_cpu_core_h']:,.0f} "
          f"core-h = {Xc['ord_cpu_core_h'] / CORES_PER_CPU_NODE:,.0f} CPU node-h")
    lab_min = [s * (44 / 29) ** 4 / 60 for s in (35.0, 107.0)]
    print(f"  optimize=True per molecule: 35-107 core-s at n29 -> {lab_min[0]:.1f}-{lab_min[1]:.1f} core-min at n44")
    # storage
    mps_gb = sum(c * S[1]['SQD_CAND'] * 32.0 * n * sample_chi(n, S[1]) ** 2 / 1e9 for n, c in SQD_SUBSET.items())
    print(f"  storage: dense MPS for all SQD states (32 n chi^2 bytes, no symmetry compression) {mps_gb / 1e3:.2f} TB; "
          f"samples {X_samples_gb(S[1]):.1f} GB; QPU bitstrings {Xc['qpu_shots'] * 16 / 1e9:.1f} GB at 16 B/shot")
    print(f"  SQD-evaluated TN states: {Xc['sqd_states']:.0f}; pretraining A100-h per full-data epoch: "
          f"{Xc['epoch_a100_h']:.3f}")

    # ---- paste-ready block: sections rounded to 10 node-h, contingency rounded, reserve = request - sum ----
    secs = [(label, round(node_h(Lc, keys) / 10.0) * 10) for label, keys in SECTIONS]
    planned = sum(v for _, v in secs)
    reserve = req - planned
    print(f"\nPaste-ready block (request {req:,} node-h = {req * 4:,} A100-h; central {tot[1]:,.0f}):")
    for l_, v in secs:
        print(f"  {l_:48s} {v:6,}")
    print(f"  {'Reserve for unbenchmarked unit costs':48s} {reserve:6,}  (request - rounded central lines)")
    print(f"  {'sum':48s} {planned + reserve:6,}")

    if md:
        print("\n---- markdown: budget table (node-hours; Low/High = MC P10/P90 per line, not additive) ----")
        print("| Component | Formula (central inputs) | Central | P10–P90 |")
        print("|---|---|--:|--:|")
        fkey = {"pre_prod": "pre_prod", "pre_sweep": "pre_sweep", "ord_labels": "ord_labels", "rl_B1": "rl_B1",
                "rl_B2": "rl_B2", "rl_B3": "rl_B3", "rl_B4": "rl_B4", "rl_heavyhex": "rl_heavyhex",
                "rl_dev": "rl_dev", "eval_tn_energy": "eval_tn_energy", "sqd_tn_prep": "sqd_tn",
                "sqd_tn_diag": "sqd_tn", "sqd_qpu_post": "sqd_qpu", "qsd": "qsd",
                "inf_timing_base": "inf_timing_base", "calib": "calib", "contingency": "contingency"}
        for label, keys in GROUPS:
            lo, _, hi = line_mc(keys)
            print(f"| {label} | {F[fkey[keys[0]]]} | {node_h(Lc, keys):,.0f} | {lo:,.0f}–{hi:,.0f} |")
        print(f"| **Total** | | **{tot[1]:,.0f}** | **{pct(tot_mc, .1):,.0f}–{pct(tot_mc, .9):,.0f}** |")
        print("\n| Section | Central | Low (P10) | High (P90) |")
        print("|---|--:|--:|--:|")
        for label, keys in SECTIONS:
            lo, _, hi = line_mc(keys)
            print(f"| {label} | {node_h(Lc, keys):,.0f} | {lo:,.0f} | {hi:,.0f} |")


def X_samples_gb(p):
    """Bytes for storing the TN samples of all SQD states (one bit per qubit, packed to bytes)."""
    return sum(c * p["SQD_SUBSET_SCALE"] * p["SQD_CAND"] * p["SHOTS"] * math.ceil(2 * n / 8)
               for n, c in SQD_SUBSET.items()) / 1e9


if __name__ == "__main__":
    main(md="--md" in sys.argv)
