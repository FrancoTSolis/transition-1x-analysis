# Task 6: n29 RL with density-matrix χ128 rewards on one GPU

*5 Oct 2026. Run: `rl_runs/followups/task6_n29_rl_dm/run/` (`log.jsonl`, `args.json`, `policy_best.pt` = step 20,
`policy_last.pt` = step 24). Driver log: `rl_runs/followups/task6_n29_rl_dm/run.log`. Energies:
`pretrain/opt_true/results/followups/task6_n29_rl_dm/`. Tables: `python3 pretrain/followups/task6_report.py`.
Every number below is measured, except the two marked as extrapolated in the last section.*

## Question

The Expanse run (`rl_runs/grpo_n29_tn`, tags `rl4n29` / `rl4n29f`) trained the n29 policy with v1 zip-up χ64
rewards on CPUs. This task asks whether a more accurate reward (the density-matrix engine at χ128, on a GPU) buys a
better policy, and what one RL step costs on one GPU.

## Answer

* **It did not buy a better policy.** At the docs' evaluation setting (v1 engine, χ256, zip margin 1.5):

  | % of CCSD correlation energy | label | start (`rl4L`) | Expanse χ64, step 15 / 30 | **this run, step 20 / 24** |
  |:--|--:|--:|--:|--:|
  | n29 val (20) | 57.2 | 63.7 | 65.7 / 65.5 | **64.7 / 62.9** |
  | n29 test (24, untouched) | 60.4 | 64.7 | 66.3 / 66.5 | **65.7 / 63.7** |

  * The best checkpoint (step 20, picked on the χ128 val) beats the start by +0.99 ± 0.30 points on test (20/24
    molecules better). That is roughly 55–60 % of the Expanse run's +1.6 to +1.8.
  * It is 0.62 ± 0.21 points below the Expanse step-15 policy on test, which is better on 17 of the 24 molecules.
  * The last checkpoint (step 24) is 0.98 ± 0.29 points below the start on test (better on only 5/24).
* **Training at χ128 was unstable.** Val at χ128 went 63.98 → 62.79 → 56.84 → 58.15 → **65.00** → 63.19
  (steps 0, 5, 10, 15, 20, 24).
  * At step 10, every val molecule was below its step-0 value (−7.1 points on average). By step 20, 18 of 20 were
    above it (+1.0).
  * On the 40 training molecules, the final policy is 1.2 points below the start (2/40 better).
  * The Expanse run never fell more than 0.3 points below its start on val. It ended +1.8 on its training set (37/40
    better).
* **The result is confounded.** This run had 6 molecules per step and 24 steps, against 12 and 30 on Expanse: 40 % of
  the Expanse run's GRPO groups. The step count was fixed in advance from the measured GPU throughput, for a
  10–12 h budget.
  * The instability looks like a noisier gradient (detail below), so batch size is the likelier cause than reward
    precision. One seed per run cannot tell the two apart.
* **χ128 changes the absolute energies but not the ranking.**
  * For the start policy on val, v1 χ64 gives 60.75 %, v1 χ256 gives 63.7 % and dm χ128 gives 63.98 %.
  * For all three policies evaluated both ways (steps 0, 20, 24), dm χ128 sits +0.28 to +0.35 points above v1 χ256 on
    average, and above it on every molecule. Gains measured at χ128 carry over to χ256 unchanged: step 20 vs 0 is
    +1.02 at χ128 and +0.95 at χ256; step 24 vs 0 is −0.80 and −0.84.
  * The policy did not learn to exploit χ128 truncation.
* **Cost.**
  * One step (48 energies) took a median 963 s (16 min) on one H200, which another user's ~96 GB job shared. That is
    20 s of GPU time per energy, or 180 energies per hour (125–250 per hour across steps).
  * The Expanse run took 262 s per step for 96 χ64 energies on one 128-core node.
  * This run took 8.0 h from driver start to the end, plus 1.3 h for the v1 χ256 evaluation. The Expanse run took
    2.5 h.

## Setup

Same as the Expanse run (`rl_runs/grpo_n29_tn/args.json`) except for the reward, the batch size, the number of steps
and the queue settings:

| | Expanse run (`grpo_n29_tn`) | this run |
|:--|:--|:--|
| reward | v1 zip-up, χ64, margin 1.5, 96 single-core CPU workers in-process | dm engine (`tn_energy.py`, `method="dm"`), χ128, complex64, through `rl_queue/task6_tn128`: 4 workers × 4 block2 threads on scai7 GPU 4 (`task6_launch.sh workers`) |
| val evals (every 5 steps) | v1 χ64 | dm χ128 |
| molecules per step × group | 12 × 8 | **6 × 8** |
| steps | 30 | **24** |
| GRPO groups seen | 360 | **144** |
| start | `grpo_large_rl4/policy_snap_for_n29.pt` (`rl4L`), 4 recycles (`--prefix-T 3`) | same |
| data | `n29_train40.txt` / `n29_val20.txt` | same |
| σ_K / σ_Z, antithetic, lr, PPO epochs, KL coef, ref update, grad clip, seed | 0.01 / 0.005, yes, 1e-5, 1, 0, 5, 1.0, 0 | same |
| driver | Expanse CPU, 8 threads | scai2 CPU, 4 threads |

* **Code.** `grpo_slot.py` got one option, `--queue-kind tn` (task kind of the queue reward). With it, the collector
  waits 1,800 s before requeuing a stale task. The default (`exact`, 300 s) is unchanged.
* **Engine pinning.** The same four worker processes served the whole run (02:05–10:20). `tn_energy.py` was last
  modified on 3 Oct, so every reward used one engine version.
* **Evaluation.** Both checkpoints were scored by `task6_eval.sh`: `policy_dump.py` on CPU, then `queue_eval.py`
  through 4 v1 χ256 workers on the same GPU after the χ128 workers were stopped (`task6_launch.sh v1workers`). The
  label, start and Expanse energies are the stored ones from `docs/rl_larger_molecules.md`.
* **Control.** Four stored `rl4n29f` v1 χ256 energies, first computed on V100/A100 GPUs, were recomputed on this H200.
  They agree to 0.21 mHa at most, about 0.05 points of the correlation energy.

## Training curve at χ128

| step | 0 | 5 | 10 | 15 | 20 | 24 |
|:--|--:|--:|--:|--:|--:|--:|
| val (20) mean % corr, dm χ128 | 63.98 | 62.79 | 56.84 | 58.15 | **65.00** | 63.19 |
| val median | 63.41 | 61.86 | 56.35 | 57.59 | 64.35 | 62.38 |
| val paired change vs step 0 (# better) | — | −1.20 (1/20) | −7.14 (0/20) | −5.84 (0/20) | +1.02 (18/20) | −0.80 (5/20) |
| train (40) mean | 62.29 | | | | | 61.10 (−1.19, 2/40 better) |

Expanse run for comparison, at v1 χ64:

| step | 0 | 5 | 10 | 15 | 20 | 25 | 30 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| val (20) mean % corr | 60.75 | 60.77 | 60.47 | **62.78** | 62.04 | 61.31 | 62.55 |
| val paired change vs step 0 | — | +0.03 | −0.28 | +2.03 | +1.29 | +0.57 | +1.81 |
| train (40) mean | 58.50 | | | | | | 60.29 (+1.79, 37/40 better) |

**Drift measure.** The training reward depends on which 6 molecules a step draws, so the raw curve is hard to read.
The batches can be replayed exactly from the driver's numpy RNG, `default_rng(0)`, which is used only for the batch
draw. Each step's mean reward minus the start policy's value on the same molecules (from the step-0 train eval at the
run's own χ) gives:

* **This run (pp):**
  * steps 1–6: −0.8, −2.5, −7.8, −8.6, −5.4, −2.1
  * steps 7–12: +0.1, −0.3, −1.1, −3.3, −5.7, −8.8
  * steps 13–18: −9.8, −10.0, −5.0, −6.6, −3.2, −2.6
  * steps 19–24: −0.6, −0.6, −0.2, +0.1, −0.5, −1.4
* **Expanse run:** between −3.4 and +1.6, and positive on steps 13–22.

Exploration noise alone costs 0.4–0.8 points: step 1 is scored before any update (−0.8 here, −0.4 on Expanse). This run fell 5–10 points below the start
twice, at steps 3–5 and 11–16, and recovered both times.

**The policy moved further each block.** The KL diagnostic measures distance from the reference policy, which resets
every 5 steps. At the end of each block it was:

* 1,005–5,782 in this run;
* 332–2,260 in the Expanse run.

That is 2.6–5.4× the Expanse value in every block. The lr is the same, and gradient clipping at norm 1 was active at
every step of both runs (norms 673–2,541 here, 514–1,646 on Expanse); only the batch is half as large. The most
likely reading is a noisier update direction from 6-molecule batches. Nothing here isolates it from the reward
change.

**No spurious rewards.** The highest single reward (71.6 %) is within 0.2 points of the best start-policy value on
the training set (71.4 %). There is no sign of unphysically low dm energies being rewarded, and no energy failed (0 of
1,152 training and 200 evaluation energies).

## χ256 (v1) results

**n29 val (20).** The best checkpoints were chosen on this set (Expanse at χ64, this run at χ128).

| candidate | mean % corr | median | min | vs label | better than label |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 57.2 | 58.7 | 40.5 | — | — |
| start: norb 15–18 RL (`rl4L`) | 63.7 | 63.2 | 57.6 | +6.5 | 19/20 |
| Expanse χ64 RL, step 15 (`rl4n29`) | 65.7 | 65.0 | 61.4 | +8.5 | 19/20 |
| Expanse χ64 RL, step 30 (`rl4n29f`) | 65.5 | 65.0 | 60.0 | +8.3 | 20/20 |
| this run, step 20 (`rl4n29dm`) | 64.7 | 64.0 | 59.6 | +7.4 | 19/20 |
| this run, step 24 (`rl4n29dmf`) | 62.9 | 62.2 | 57.8 | +5.7 | 18/20 |

**n29 test (24, untouched).**

| candidate | mean % corr | median | min | vs label | better than label |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 60.4 | 61.1 | 48.6 | — | — |
| start (`rl4L`) | 64.7 | 64.3 | 56.5 | +4.3 | 22/24 |
| Expanse step 15 (`rl4n29`) | 66.3 | 66.0 | 58.7 | +5.9 | 24/24 |
| Expanse step 30 (`rl4n29f`) | 66.5 | 66.7 | 59.5 | +6.1 | 24/24 |
| this run, step 20 (`rl4n29dm`) | 65.7 | 65.9 | 57.5 | +5.2 | 23/24 |
| this run, step 24 (`rl4n29dmf`) | 63.7 | 64.2 | 54.5 | +3.3 | 20/24 |

**Test by formula** (mean % corr, 6 molecules each). C3H5NO2 and C5H5N appear in no RL set.

| formula | label | start | Expanse 15 | Expanse 30 | this run 20 | this run 24 |
|:--|--:|--:|--:|--:|--:|--:|
| C3H5N3 | 56.9 | 62.4 | 64.2 | 64.3 | 63.5 | 61.5 |
| C3H5NO2 | 62.5 | 65.6 | 66.4 | 66.6 | 65.2 | 63.6 |
| C4H5NO | 58.4 | 63.8 | 65.9 | 65.7 | 64.8 | 62.7 |
| C5H5N | 63.9 | 66.9 | 68.6 | 69.3 | 69.2 | 66.9 |

**Paired differences** (percentage points, mean ± s.e.m., # molecules where the first is better):

| | val (20) | test (24) |
|:--|:--|:--|
| this run step 20 − start | +0.95 ± 0.22 (18) | +0.99 ± 0.30 (20) |
| this run step 24 − start | −0.84 ± 0.23 (3) | −0.98 ± 0.29 (5) |
| this run step 20 − Expanse step 15 | −1.05 ± 0.21 (2) | −0.62 ± 0.21 (7) |
| this run step 20 − Expanse step 30 | −0.80 ± 0.13 (1) | −0.83 ± 0.18 (4) |
| this run step 24 − Expanse step 30 | −2.59 ± 0.18 (0) | −2.80 ± 0.27 (0) |

* **By formula.** Step 20 matches the Expanse policies on C5H5N (69.2 against 68.6 / 69.3). It trails them by
  0.7–1.4 points on the other three formulas, including C3H5NO2, the other formula outside the RL sets. Neither
  out-of-set formula shows a consistent difference in how the two runs generalize.
* **Use of the val set.** Picking step 20 used the val set at χ128, so its val numbers are biased upward (as are those
  of `rl4n29`). The test numbers are not.

## Cost

| | this run (dm χ128) | Expanse run (v1 χ64) |
|:--|:--|:--|
| hardware | one H200 (scai7 GPU 4) shared with another user's 96.4–96.7 GB job; 16 scai7 cores for the block2 ⟨H⟩ threads; driver on 4 scai2 cores | one 128-core Expanse node (96 single-core workers + 8-thread driver) |
| energies per step | 48 | 96 |
| step wall clock, median (range) | 963 s (691–1,381) | 262 s (256–268) |
| per energy | 20.0 s of the GPU's wall clock (180 / h) | 2.7 s of the node's wall clock = 260 core-s |
| val eval (20 energies) | 541–571 s (27–29 s per energy) | 137–147 s |
| train eval (40) | 573–1,144 s (14–29 s per energy) | 150–152 s |
| RL steps, total | 6.60 h (24 steps) | 2.18 h (30 steps) |
| whole run | 7.85 h from first val eval to end; 8.0 h from driver start | 2.50 h |

* **GPU footprint** (nvidia-smi, 5-s samples by `gpu_monitor.py`). Both phases stayed far below the 40 GB cap.
  * The four dm χ128 workers used 8.8 GB (median) and 13.5 GB at most, with at most 5.6 GB per process.
  * The four v1 χ256 evaluation workers used 6.9 GB (median) and 8.0 GB at most, with at most 3.2 GB per process.
  * Mean GPU utilization was 93 % during RL and 72 % during evaluation. These figures include the other user's job,
    which kept the GPU near 100 % on its own much of the time.
* **v1 χ256 evaluation.** 92 energies took 1.29 h on 4 workers: 2 checkpoints × (20 + 24), plus the 4-energy control.
  Each energy took a median of 206 s (142–263 s) of worker time. That is ~72 energies per hour on this GPU.
* **The per-energy cost varied with the shared GPU.** Training steps ran 125–250 energies per hour. The two 40-energy
  train evals, on the same molecules, took 1,144 s and 573 s. The workers do not log per-task phases, so the spread
  cannot be assigned. These numbers are therefore not a clean benchmark.
  * The median (180 per hour) is above the 126 per hour measured in the task-4 queue test on a shared H200.
  * The val evals (27–29 s per energy, ~130 per hour) match task 4.

## What the χ128 reward bought

* **More accurate absolute energies.** dm χ128 scores policies about 3.2 points closer to convergence than the v1 χ64
  reward, and 0.3 points closer than the docs' own v1 χ256 evaluation.
  * This matters for reporting energies.
  * RL uses only the ranking within a group of 8, and χ64 already got that right: median Spearman 0.983 in
    `rl_larger_molecules.md`. The χ64 run's gains survived at χ256 unchanged.
* **No better policy at this budget.**
  * With 40 % of the Expanse run's groups, the χ128 run's best checkpoint gained about 1 point on the untouched test,
    against 1.6–1.8 for the χ64 run.
  * Its last checkpoint lost 1 point.
  * Per GRPO group (one molecule × 8 samples), a χ128 step cost 2.7 minutes of a shared H200 (963 s / 6 groups).
    An Expanse χ64 step cost 22 s of a full node (262 s / 12 groups).
* **Not separated: batch size and reward.** The instability and the smaller gain fit the halved batch, but one seed
  per run cannot rule out a reward effect. Two ways to separate them:
  * a v1 χ64 run with 6 molecules per step, about 2 h on one Expanse node (**extrapolated**: each step is still one
    round of the 96 workers, 262 s, × 24 steps, plus evals);
  * a dm χ128 run with the Expanse recipe (12 × 8, 30 steps). On this shared H200 that would take about 16 h
    (**extrapolated**: 2 × 963 s per step × 30 steps, assuming the reward time scales linearly with the number of
    energies), or about 8 h on two such GPUs.
* **For the NERSC RL line.** At 29 orbitals nothing here argues for χ128 over χ64 rewards. Group-ranking fidelity at
  33–44 orbitals, where χ64's error against χ256 grows to 0.42–1.15 points (task 5), is still untested.

## Files

* **Run.** `rl_runs/followups/task6_n29_rl_dm/`:
  * `run/` (`args.json`, `log.jsonl`, `policy_best.pt` = step 20, `policy_last.pt` = step 24);
  * `run.log`;
  * `eval.log`, `eval_{control,val,test}.log`;
  * `tasks_n29{val,test}_rl4n29dm.pkl` (the energy tasks of both checkpoints);
  * `gpu4_monitor_partB.jsonl` (GPU footprint);
  * `worker_pids.txt`, `v1_worker_pids.txt`;
  * `control_val4.txt`.
* **Energies.** `pretrain/opt_true/results/followups/task6_n29_rl_dm/`:
  * `energy_n29val_rl4n29dm_v1chi256.json` and `energy_n29test_rl4n29dm_v1chi256.json` (tags `rl4n29dm` = step 20 and
    `rl4n29dmf` = step 24);
  * `energy_n29val_control_v1chi256.json`;
  * `summary_pairs.json` (paired differences, χ128 vs χ256 offsets, control).
* **Code.**
  * `pretrain/followups/task6_launch.sh` (χ128 workers, driver, v1 workers);
  * `task6_eval.sh`;
  * `task6_report.py` (every table above, plus the batch replay and cost and footprint summaries);
  * `gpu_monitor.py`;
  * `pretrain/rl/grpo_slot.py` (`--queue-kind`).
* **Processes.** All stopped:
  * the driver exited normally at 10:19;
  * the χ128 workers were stopped at 10:20;
  * the v1 workers and the GPU monitor were stopped at 11:39.

  scai7 GPU 4 now holds no fts process.
