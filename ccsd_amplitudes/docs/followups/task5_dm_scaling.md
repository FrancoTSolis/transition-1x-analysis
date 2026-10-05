# Task 5: the density-matrix TN engine at 33–44 orbitals (accuracy, cost, and the NERSC RL line)

*5 Oct 2026. One H200 on scai7 (GPU 4), which another user's inference server shared the whole time. Results:
`pretrain/opt_true/results/followups/task5_dm_scaling/` (`bench.jsonl`, `queue_tput.json`, `summary.json`). Logs:
`rl_runs/followups/task5_dm_scaling/`. Code: `pretrain/followups/task5_*.py`, `gpu_monitor.py`. Every table below is printed by
`python3 pretrain/followups/task5_analyze.py`.*

## Question

Nobody had run the engine above 29 orbitals. How accurate are its energies (`pretrain/rl/tn_energy.py`,
`method="dm"`) at 33, 37 and 44 orbitals for χ = 64, 128 and 256? What do they cost? Do the cost model's assumptions
hold (`docs/nersc_2027_cost_model.py`: MPS-build exponent 2.6 in n, ⟨H⟩ exponent 4.0, χ exponent 1.8, 4 workers
sharing one GPU)? Is the RL line of the 10,000-node-hour request (2,800 node-hours) still consistent?

## Answer

* **Accuracy holds up at 44 orbitals.** At 33–44 orbitals:
  * χ128 sits 0.28–1.81 mHa above χ256 (median 1.06), which is 0.06–0.41 points of the CCSD correlation energy;
  * χ64 sits 1.9–5.2 mHa above χ256 (0.42–1.15 points);
  * at 29 orbitals the same gaps are 0.62–0.69 and 2.5 mHa.

  On the one χ512 run (C7H16 TS, 44 orbitals), χ256 is 0.68 mHa (0.13 points) above χ512. A linear extrapolation of E
  in the discarded weight through χ128/χ256 had predicted that χ512 energy to within 0.04 mHa. The energy error per
  unit discarded weight is 1.2–2.4 Ha, inside the 0.9–6.2 Ha calibrated against exact energies at norb 15–18.
  `discarded_sum` therefore stays a usable error bar: the estimated χ256 residual is 0.1–1.1 mHa at 33–44 orbitals.
  TS geometries are the hardest cases.
* **Cost grows more slowly with size than the model assumes.** One worker, quiet GPU, χ128, 29 → 44 orbitals:
  * MPS build ∝ n^2.07±0.13 (21.6 → 50 s; model 2.6);
  * block2 ⟨H⟩ ∝ n^3.98±0.04 (4.0 → 20.9 s on 4 Zen 5 threads; model 4.0);
  * MPO build ∝ n^5.72±0.02 (6.3 → 68 s; model 4, divided by 8).
* **The 4 workers on one GPU do not run their GPU work concurrently.** This GPU time-slices CUDA contexts. With 2–3
  other builds running, a build took 1.9–4.4× longer, while ⟨H⟩ (CPU) was unchanged. Through the real queue (4
  workers × 4 threads, χ128, the RL configuration) the throughput was:

  | orbitals | energies per GPU-hour | ratio to 29 orbitals | single worker |
  |:--|--:|--:|:--|
  | 29 | 126 (task 4) | 1.00 | 1.1× slower |
  | 37 | 107 | 1.18 | – |
  | 44 | 53 | 2.36 | 2.0× slower |

  The model's ratio from 29 to 44 orbitals is 3.77. A GPU-bound model, `max(build, (build + CPU phases)/4)` with the
  clean single-worker phase times, predicts 2.37.
* **The RL line is consistent.**
  * With the measured growth and the model's own 29-orbital anchor (34.8 A100-s per energy), the RL lines come to
    2,390–2,920 node-hours (2,800 budgeted).
  * With the anchor set to the measured 4-worker throughputs (H200 28.7 s, TITAN Xp 48.8 s per energy, task 4), they
    come to 1,980–4,070.
  * The size scaling does not raise the request.
  * The open risks are the per-GPU throughput at 29 orbitals on an A100 (only bracketed so far) and two engineering
    fixes (below) that decide whether production sits at the low or the high end.
* **Memory (nvidia-smi, per process).**

  | χ | footprint |
  |:--|:--|
  | 64 | 1.2–1.5 GB |
  | 128 | 2.4–6.1 GB |
  | 256 | 5.8–8.1 GB |
  | 512 (44 orbitals) | 15.5 GB |

  * The real footprint is torch's reserved pool plus 0.7–0.9 GB.
  * Torch's allocated peak is 1.3–5× smaller than the reserved pool, worst at χ512 because of fragmentation.
  * Four χ256 workers fit a 40 GB A100 (≤ 33 GB).

## Method

* **Molecules.** Three per size: one R, one P and one TS from three different reactions (seeded draw), plus two
  29-orbital anchors that appear in no RL list. The molecule lists are in `names.txt`. Hamiltonians were built with
  `pretrain/rl/hamiltonian.py`; E_HF and E_CCSD reproduce `rhf_dataset` exactly. Frames exist for every shape.

| orbitals | formula (split in the cost model) | (nocc, nvirt) | molecules | E_corr(CCSD), mHa |
|--:|:--|:--|:--|:--|
| 29 | C3H5N3, C4H5NO | (16, 13) | rxn2003_P, rxn3873_P | −431, −399 |
| 33 | C5H5NO (test) | (18, 15) | rxn7064_R, rxn7940_P, rxn5687_TS | −479, −451, −494 |
| 37 | C4H9NO2 (val) | (21, 16) | rxn5799_R, rxn7572_P, rxn6579_TS | −439, −435, −438 |
| 44 | C7H16 (test) | (22, 22) | rxn8743_R, rxn8744_P, rxn8742_TS | −488, −475, −526 |

* **Parameters.** The `rl4n29f` policy (`rl_runs/grpo_n29_tn/policy_last.pt`), dumped with `policy_dump.py` on CPU.
  It scores 66–79 % of the CCSD correlation energy here. For a same-input check against the TITAN Xp runs of
  `tests/results/tn_v2`, the `net` parameters of C3H5N3_rxn2003_P were also run.
* **Engine.** `LUCJEnergyTN` with:
  * method `dm`, complex64, cutoff 1e-9, mode_tol 1e-8;
  * Boys split-localized basis, Fiedler order;
  * block2 ⟨H⟩ on 4 OpenMP threads, stack 2 GB;
  * one engine per molecule per process, so the MPO is built once per molecule and its build time is recorded on that
    molecule's first job.

  Driver: `pretrain/followups/task5_dm_bench.py`, launched by `task5_launch.sh`. Every process had all five thread
  variables set to 4 and a torch memory cap (9 GB in pass 1, 30 GB later).
* **Hardware.** scai7 GPU 4 (H200, 141 GB) and an AMD EPYC 9535 host (Zen 5). The GPU was shared with another user's
  `sglang` server, which held 97 GB and was ≥ 95 % busy for most of the night. That load is the main noise source in the
  GPU timings.
* **Passes.** 66 energies in total, none failed; repeated jobs reproduced E to ≤ 1e-6 Ha.
  * Pass 1: 4 concurrent processes (the cost model's 4 workers per GPU), 33 molecule/χ jobs plus the TITAN-check jobs.
  * Pass 2: one process alone, 21 jobs, including χ512.
  * Pass 3: automatic retries of busy runs; stopped early because the GPU stayed busy.
* **Clean timings.** A build counts as clean if no other build of mine overlapped it and the whole GPU was ≥ 95 % busy
  for at most 25 % of it. One worker alone keeps the GPU 7–50 % busy (mean utilization during its build), so higher
  values mean the foreign job was active. CPU phases (⟨H⟩, MPO) are unaffected by GPU sharing (ratio 0.99–1.05), so their medians use all runs.
* **Memory.** `pretrain/followups/gpu_monitor.py` logged the per-process `nvidia-smi` footprint once per second
  (`gpu4_monitor_partA.jsonl`).
* **Queue throughput.** 4 standard workers: `reward_queue worker --kind tn --chi 128 --tn-impl current
  --block2-threads 4`, THREADS=4, on `rl_queue/task6_tn128`; the same workers then served task 6. The script
  (`task5_queue_tput.py`) submitted 8 tasks at 44 orbitals and then 8 at 37 orbitals, in random order. As in a GRPO
  step, most tasks rebuild their MPO.

## Accuracy

% of the CCSD correlation energy, ΔE against the largest χ (mHa, positive = above), and the summed discarded weight
of all truncations. "χtop − extrap" is the remaining error of the largest χ from a linear fit of E against the
discarded weight through the two largest χ.

| molecule | norb | % χ64 | % χ128 | % χ256 | % χ512 | ΔE χ64 | ΔE χ128 | ΔE χ256 | disc χ128 | disc χ256 | χtop − extrap | slope (Ha) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| C3H5N3_rxn2003_P | 29 | 68.02 | 68.45 | 68.60 | | 2.47 | 0.62 | 0 | 4.7e-4 | 1.3e-4 | 0.22 | 1.79 |
| C4H5NO_rxn3873_P | 29 | 68.46 | 68.91 | 69.09 | | 2.49 | 0.69 | 0 | 4.9e-4 | 1.3e-4 | 0.24 | 1.92 |
| C5H5NO_rxn5687_TS | 33 | 65.68 | 66.44 | 66.73 | | 5.19 | 1.44 | 0 | 1.2e-3 | 3.7e-4 | 0.66 | 1.79 |
| C5H5NO_rxn7064_R | 33 | 70.87 | 71.30 | 71.44 | | 2.72 | 0.66 | 0 | 7.7e-4 | 2.8e-4 | 0.39 | 1.37 |
| C5H5NO_rxn7940_P | 33 | 71.34 | 71.68 | 71.76 | | 1.92 | 0.39 | 0 | 4.6e-4 | 1.4e-4 | 0.17 | 1.23 |
| C4H9NO2_rxn5799_R | 37 | 72.05 | 72.65 | 72.93 | | 3.87 | 1.23 | 0 | 9.4e-4 | 3.2e-4 | 0.62 | 1.97 |
| C4H9NO2_rxn6579_TS | 37 | 67.83 | 68.56 | 68.97 | | 5.02 | 1.81 | 0 | 1.2e-3 | 4.3e-4 | 1.04 | 2.41 |
| C4H9NO2_rxn7572_P | 37 | 70.47 | 70.86 | 70.92 | | 1.95 | 0.28 | 0 | 3.1e-4 | 9.0e-5 | 0.11 | 1.24 |
| C7H16_rxn8742_TS | 44 | 77.52 | 78.03 | 78.29 | **78.42** | 4.72 | 2.03 | 0.68 | 1.5e-3 | 5.1e-4 | 0.34 (χ512) | 2.01 |
| C7H16_rxn8743_R | 44 | 78.86 | 79.22 | 79.35 | | 2.39 | 0.61 | 0 | 6.1e-4 | 1.9e-4 | 0.27 | 1.44 |
| C7H16_rxn8744_P | 44 | 78.79 | 79.49 | 79.72 | | 4.42 | 1.06 | 0 | 1.2e-3 | 5.2e-4 | 0.78 | 1.51 |

* For C7H16_rxn8742_TS the ΔE values are relative to χ512 (discarded weight 1.7e-4). Its χ128 → χ256 gap is
  1.35 mHa.
* The same-input check reproduces the TITAN Xp energies of C3H5N3_rxn2003_P (`net`) to 0.004 / 0.023 / 0.055 mHa at
  χ64 / 128 / 256.
* All errors are positive, as expected from truncation. They follow the discarded weight, which is largest for TS
  geometries, more than they follow the size.
* For RL, χ128 at 33–44 orbitals is about as accurate as χ128 at 29, the setting validated for rewards. Group-ranking
  fidelity at these sizes was not measured.

## Cost per energy

**Clean single-worker times** (s; medians over molecules; nvidia-smi footprint in MiB):

| norb | χ | molecules | MPS build | (env + sweep) | ⟨H⟩, 4 threads | MPO build | footprint |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 29 | 128 | 1 | 21.6 | 6.7 + 14.8 | 4.0 | 6.3 | 3,898 |
| 33 | 64 | 1 | 12.8 | 2.3 + 10.4 | 3.4 | 13.0 | 1,182 |
| 33 | 128 | 2 | 26.4 | 8.7 + 17.6 | 6.6 | 13.0 | 2,701 |
| 33 | 256 | 1 | 128.7 | 54.7 + 73.8 | 17.4 | 13.0 | 5,830 |
| 37 | 64 | 1 | 14.0 | 2.7 + 11.1 | 5.3 | 25.2 | 1,290 |
| 37 | 128 | 2 | 32.2 | 10.9 + 21.1 | 10.3 | 25.2 | 5,495 |
| 37 | 256 | 1 | 132.2 | 55.6 + 76.4 | 26.8 | 25.2 | 7,406 |
| 44 | 128 | 2 | 50.0 | 15.1 + 34.7 | 20.9 | 67.7 | 5,225 |
| 44 | 512 | 1 | 696 (GPU busy) | 360 + 335 | 136.6 | 67.7 | 15,530 |

Missing rows had no quiet run.

**Exponents.**

| quantity | measured | cost model |
|:--|:--|:--|
| MPS build vs n at χ128 (7 clean runs, 29–44 orbitals) | 2.07 ± 0.13 | A_STATE 2.6 (range 2.0–3.2) |
| ⟨H⟩ vs n | 3.98 ± 0.04 (3.8–4.0 at other χ) | A_EXPECT 4.0 |
| MPO build vs n (20 builds) | 5.72 ± 0.02 | A_MPO 4, amortized over 8 |
| MPS build vs χ, 64 → 128 (clean, 33–37 orbitals) | 1.2 | 1.2 |
| MPS build vs χ, 128 → 256 | 2.0–2.3 | 1.8 |
| ⟨H⟩ vs χ, 64 → 128 / 128 → 256 / 256 → 512 | 0.97 / 1.37 / 1.38 | – |

**Sharing one GPU.** Ratios of pass-1 build time (4 processes) to the clean single-worker time for the same jobs:
1.9–4.4 at χ64–256 with 2–3 concurrent builds. ⟨H⟩ is unchanged (ratio 0.99–1.03). Same inputs on the TITAN Xp
(loaded Broadwell host): build 94 s at χ128 and ⟨H⟩ 122 s on 1 thread. Here, for the same molecule with the rl4n29f
parameters, the clean build took 21.6 s and ⟨H⟩ 4.0 s on 4 threads.

**4-worker throughput through the queue** (the RL configuration; most tasks rebuild their MPO):

| norb | tasks | wall (s) | s per energy | energies / GPU-h | ratio to 29 | GPU-bound model (s) | one worker, MPO every task (s) |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 29 | 40 | 1,147 | 28.7 | 126 | 1.00 | 20.4 | 30.7 |
| 37 | 8 | 270 | 33.8 | 107 | 1.18 | 33.8 | 69.5 |
| 44 | 8 | 540 | 67.5 | 53 | 2.36 | 48.4 | 137 |

* The 29-orbital row is the task-4 test, run on the same GPU at 00:10, before this task. Model columns use the fitted
  clean phase laws.
* One of the 16 tasks was computed twice (an orphan result). This is the known requeue race of `reward_queue.py`; it is
  harmless because the first result is kept.

## The NERSC cost model and the RL line

* **The model's structure.** Its per-energy cost is `[build + ⟨H⟩ + MPO/8] × 1.2 / 4`. Here the build is GPU-bound and
  serialized across the 4 processes, so it does not divide by 4.
* **At 29 orbitals the anchor still roughly holds, for compensating reasons.** The model's build time (75.5 A100-s,
  half the time on a loaded TITAN Xp host) is ~3.5× longer than the clean H200 + Zen 5 build (21.6 s). Its division by
  4 workers is about 4× too optimistic on this GPU. The two errors largely cancel. The measured 4-worker cost per energy
  brackets the model's 34.8 A100-s: 28.7 s on the shared H200 and 48.8 s on a TITAN Xp (task 4).
* **What task 5 adds is the growth with size.** The table applies the measured growth to the model's buckets and keeps
  every other model input:

| growth with n (29 → 44 ratio) | anchor 28.7 s (H200) | 34.8 s (model) | 48.8 s (TITAN Xp) |
|:--|--:|--:|--:|
| model (3.77) | 2,320 | **2,800** | 3,910 |
| measured, GPU-bound (2.37; reproduces the queue data) | 1,980 | **2,390** | 3,330 |
| measured phases, all parallel, MPO/8 (3.08) | 2,140 | 2,580 | 3,600 |
| measured phases, all parallel, MPO every task (4.47) | 2,420 | 2,920 | 4,070 |

* All numbers are RL node-hours: the B1 exact line plus the TN lines (B2–B4, heavy-hex, development). The model row
  scales the 2,800 by the anchor.
* B3 (30–37 orbitals) and B4 (38–44), 74 % of the RL line, drop from 1,004 / 1,075 to 858 / 781 node-hours under the
  GPU-bound growth at the model's anchor.

**Verdict.** The size scaling beyond 29 orbitals does not raise the RL line. At the model's anchor the line comes to
2,390–2,920 node-hours, against 2,800 budgeted. The RL line moves with the A100 per-energy cost at 29 orbitals more
than with anything measured here:
* if an A100 behaves like this H200, the line is ~2,000–2,400;
* if it behaves like a TITAN Xp, it is ~3,300–4,100. That would use up most of the 17.5 % contingency on its own.

Two fixes decide where in that range production lands, and both are cheap:

1. **Let the 4 workers share the GPU.** Run them under the CUDA MPS daemon, which NERSC allows in user jobs, or build
   the 8 states of a GRPO group in one process. A single worker keeps the H200 busy only 7–50 % of the time, so either should
   recover a large part of the 4× the model assumes.
2. **Reuse the MPO.** The queue workers rebuild the block2 MPO for almost every task, because tasks of different
   molecules arrive in random order and `--max-items 1` is used. The MPO build grows as n^5.7: at 44 orbitals it takes
   68 s, more than the 50 s MPS build. In today's GPU-bound regime this is hidden behind the other workers' builds. Once
   fix 1 lands it would dominate, unless tasks are routed by molecule (the model assumes one MPO per 8-sample group).

## Caveats

* **Measured on an H200 shared with a foreign job, not on an A100.** Clean GPU timings exist only where the foreign
  job was idle: 7 χ128 runs, and only 2 sizes at χ64 and χ256. The 44-orbital χ256 and χ512 timings, the
  29-orbital χ64 and χ256 timings, and every pass-1 number are inflated. The exponents are ratios on one machine; the
  A100 anchor is not measured.
* **A busy proxy.** The "quiet" test uses whole-GPU utilization. nvidia-smi does not report the other user's
  per-process load, and my own small kernels can also read as busy (29 orbitals, χ64).
* **One molecule type per size.** Each size has one formula. The 44-orbital C7H16 (an alkane, 78–80 %) is easier per
  orbital than the N/O molecules, so the 44-orbital accuracy may be optimistic for other formulas at that size.
  Accuracy per size rests on 3 states of one policy, and χ512 on a single state.
* **The 4-worker queue test is short.** 8 tasks per size, one run each, during the foreign job.
* **The CPU side differs from Perlmutter.** ⟨H⟩ and the MPO ran on Zen 5 cores at up to 4.3 GHz; Perlmutter's EPYC
  7763 is slower per core. The cost model's ⟨H⟩ anchor (35 s at 4 threads, from scai2) is far above the 4.0 s measured
  here.

## What to run next

1. On a Perlmutter A100 (or the scai3 A100 once reachable), with no other job on the GPU:
   * the clean single-worker build at 29 and 44 orbitals, χ128;
   * 4-worker queue throughput with and without the CUDA MPS daemon.

   That pins the anchor that dominates the RL line.
2. Molecule-affine task claiming, or larger per-worker engine caches, in `reward_queue.py`, so each MPO is built once
   per GRPO group (the cost-model assumption).
3. Group-ranking fidelity at χ128 for 37–44-orbital perturbation groups, against χ256/χ512. Only the absolute accuracy
   was measured here.
4. A χ512/χ1024 memory test with `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` for the sampling MPS of the
   TN-SQD tier. At χ512 the reserved pool was 5× the allocated peak.
