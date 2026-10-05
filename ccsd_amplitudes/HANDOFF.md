# Handoff: ML-LUCJ (ccsd_amplitudes/), state as of 5 Oct 2026

Read this first, then `docs/rl_larger_molecules.md` (latest results), `docs/chem_frame_pretraining.md` (how the
pretraining works) and `docs/optimize_true_learnability_study.md` (why the old label regression failed).

## 1. Where things stand

* **What exists.** A network predicts LUCJ parameters (2 layers, square connectivity, spin-balanced; U and Z/J; the
  final orbital rotation comes from t1) from CCSD amplitudes of Transition1x molecules (STO-3G, frozen core,
  15–44 orbitals). The model is a chemistry-frame slot model pretrained on the optimize=True compressed-DF objective
  (`runs_ot/slotall_T4/best.pt`, 4 recycles), then fine-tuned with GRPO on the LUCJ variational energy.
* **Headline.** One policy (`rl_runs/grpo_n29_tn/policy_last.pt`, tag `rl4n29f`) beats the per-molecule
  optimize=True labels at every size tested, in % of the CCSD correlation energy:

  | set | energies | labels | policy |
  |:--|:--|--:|--:|
  | norb 15–16 (9 molecules) | exact | 62.0 | 70.3 |
  | norb 17–18 (22) | exact | 61.7 | 68.6 |
  | norb 19 (12, untouched) | exact | 66.7 | 68.9 |
  | n29 test (24, untouched, 4 formulas) | MPS χ256 | 60.4 | 66.5 |

* **New (5 Oct, zero-cost analysis).** On 10 complete R/P/TS reactions at norb 17–19, all LUCJ energies overestimate
  barriers. Against CCSD(T), the barrier MAE is 14.5 kcal/mol for the policy, 17.0 for the labels, 20.8 for
  Hartree–Fock and 4.6 for CCSD (`pretrain/rl/reaction_energies.py`,
  `pretrain/opt_true/results/reaction_energies_n17_19.json`). 2-layer LUCJ recovers less correlation at TS geometries.
* **NERSC.** The 2027 ERCAP GPU budget note for Brandon is `docs/nersc_2027_gpu_budget.md`: 10,000 node-hours, built
  by `docs/nersc_2027_cost_model.py` (run it with plain python3; every input is commented with its source). The
  deadline appears to be 5 Oct 2026, 5 pm PDT.

## 2. Code map

| path | what |
|:--|:--|
| `pretrain/opt_true/` | frames, slot models, `train_slot_all.py`, `eval_slot*.py`, `dftorch.py` (compressed-DF objective) |
| `pretrain/rl/grpo_slot.py` | GRPO driver. `--reward exact\|tn\|queue` (queue = remote GPU workers); `--init-policy`, `--prefix-T 3` (frozen recycles) |
| `pretrain/rl/policy_dump.py` | candidates → energy-task pickle. Tags: `init:`, `init4:` (pretrained), policy checkpoints; `--labels-dir` adds the label |
| `pretrain/rl/gpu_energy.py` | exact LUCJ energy on one GPU, complex64 to ~1e-7 Ha. Up to norb 18 on 40 GB, norb 19 on ≥ 60 GB. `state()` returns the CI matrix (for sampling) |
| `pretrain/rl/tn_energy.py` | MPS energy, density-matrix engine (`method="dm"` on CUDA, `"zipup"` on CPU); verified |
| `pretrain/rl/tn_energy_v1.py` | frozen zip-up engine. **Every n29 number in the docs used it** (`--tn-impl v1`) |
| `pretrain/rl/reward_queue.py`, `launch_workers.sh` | shared-filesystem energy queue; GPU guard (≤ 2 user GPUs per machine on scai3–7); `worker --kind exact\|tn --chi --tn-impl --zip-margin --block2-threads`; env `THREADS`, `WORKER_ARGS` |
| `pretrain/rl/queue_eval.py` | evaluate a task pickle through a queue (`--kind tn`, `--only-names`, `--skip-keys`) |
| `pretrain/rl/report_tables.py` | prints every table of `docs/rl_larger_molecules.md` |
| `pretrain/rl/hamiltonian.py`, `generate_compressed_targets.py` | active-space Hamiltonians (CPU, seconds); optimize=True labels (`--maxiter`, default 500) |
| `pretrain/rl/ccsd_t_refs.py`, `reaction_energies.py` | CCSD(T) references (cached JSON); reaction energies and barriers |
| `expanse/` | `expanse.py` (TOTP + ControlMaster helper: `run`/`put`/`get`), `sync.sh` (+ `sync_extra.txt`), slurm scripts (`grpo_n29_tn.slurm`, …) |

**Evaluation sets.** `pretrain/rl/small_val.txt` (9, norb 15–16), `large_val.txt` (22, norb 17–18),
`gauge_study/names_norb19.txt` (12), `pretrain/rl/n29_val20.txt`, `n29_test24.txt` (untouched; C3H5NO2 and C5H5N are
in no RL set).

**Candidate tags.** `label`, `pre1`/`pre4` (pretrained, 1 or 4 recycles), `rl1s`/`rl4s` (small-set RL),
`rl4L` (`grpo_large_rl4` step 25), `rl4Lf` (step 50), `rl4n29`/`rl4n29f` (`grpo_n29_tn` best/last).

**Results.** JSON in `pretrain/opt_true/results/energy_*.json`; RL logs in `rl_runs/*/log.jsonl`. Checkpoints (`*.pt`),
labels (`rhf_targets_compressed*`), Hamiltonians and queues are gitignored.

## 3. Compute rules (from Fang; keep)

* **scai2** (local, 8× TITAN Xp 12 GB, 48 cores): free to use, except **never GPUs 0/1**.
* **scai1** (8× V100 16 GB, 80 cores): free, but skip busy GPUs. On 5 Oct it was saturated by other users (load ~165,
  all GPUs busy).
* **scai3–scai7**: at most **2 GPUs per machine for user fts, counting all of Fang's jobs**. scai4 and scai6 are
  usually already at 2 (his gpu4pyscf jobs), so they are off-limits. scai7 H200s are shared with other users'
  70–130 GB jobs: cap memory (`--max-mem-gb`). scai3 ssh hung on 5 Oct (ping OK).
* **CPU-heavy work** only on scai1/scai2 (or Expanse).
* **Expanse** (account cla361, shared): keep jobs economical. Use `expanse/expanse.py`. Never print or copy the TOTP
  seed, never touch the remote `fermi_arc/` directory or the conda env `negf`, and never touch the original
  batch_screen project.
* Never kill processes you did not start.

## 4. Pitfalls already hit

* **Lazy imports.** TN workers import the engine on their first task. Swapping `tn_energy.py` mid-run silently mixes
  engine versions. Pin `--tn-impl` for any result set.
* **Heartbeats.** They come from a child process (block2 holds the GIL; a heartbeat thread starved). TN collectors use
  `stale_s=1800`. Old orphan results in `rl_queue/*/done` are harmless.
* **Backgrounding.** `cmd | tail && nohup x &` backgrounds the whole `&&` list. This once submitted an evaluation twice.
* **GPU memory.** An exact norb-18 worker next to TN workers on a 40 GB A100 ran out of memory. Size workers by the
  nvidia-smi footprint, which is 2–3× torch's allocated memory.
* **Absolute values.** MPS χ256 (v1, zip margin 1.5) is 0.4–1.5 points below converged at n29. Comparisons between
  candidates were stable to < 0.1 point across χ128 / χ256 / χ256 at margin 3.
* **Labels.** They stop at the 500-iteration L-BFGS cap; almost none converge (task 3 below).

## 5. Open tasks (Fang, 5 Oct; ordered by value)

| # | task | status / plan |
|:--|:--|:--|
| 1 | Noiseless SQD comparison with Lin et al.'s protocol: exact GPU state vectors at 15–19 orbitals; 10⁵ samples, 10 batches × 4,000, her SQD settings; candidates truncated CCSD init / label / pretrained / NN+RL; references CCSD(T), FCI at 15–16 | launched 5 Oct (workflow agent) |
| 2 | Per-molecule direct optimization, 500 NOMAD/SPSA evaluations, 15–16 orbitals. Objective: LUCJ energy (RL upper bound) and QSCI energy (her TN-optimization baseline) | launched 5 Oct |
| 3 | optimize=True labels run to convergence (~5,000 iterations) on the test sets, then energies (exact ≤ 19; TN v1 χ256 at 29) | launched 5 Oct |
| 4 | dm TN engine timing, 29 orbitals, χ128, 4 workers × 4 threads per GPU (the NERSC RL line) | running 5 Oct 00:10 on a TITAN Xp (clean host) and an H200 (shared); scai3 A100 unreachable |
| 5 | dm engine accuracy and timing at 33 / 37 / 44 orbitals, χ 64 / 128 / 256 | launched 5 Oct |
| 6 | n29 RL redone with dm χ128 rewards on GPUs | launched 5 Oct (after task 5) |
| 7 | more seeds of the norb 15–18 exact-reward RL | launched 5 Oct |
| 8 | reaction energies and barriers from existing energies | **done** (§1) |
| – | MPS sampler for TN-sampled SQD at 29 orbitals (~1 day of code) | launched 5 Oct (development) |
| – | GPU SQD solver (SBD) build and benchmark | not started |

Results of the follow-up tasks go to `docs/followups_2026-10.md` and `pretrain/opt_true/results/followups/`.
