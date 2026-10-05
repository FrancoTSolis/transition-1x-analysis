# Task 3: optimize=True labels run to convergence

5 Oct 2026. Fits by the first task-3 agent (00:20–02:05); n29 TN energies finished at 04:02; analysis rerun and report
written by the resuming agent (09:45). Every number below is measured; nothing is extrapolated.

## Answer

* **Running the labels to convergence does not change the label baseline.** The labels behind every table in
  `docs/rl_larger_molecules.md` stop at the 500-iteration L-BFGS cap; only 12 of 87 had converged. Run to scipy's
  L-BFGS-B convergence (all 87 converge, after 401–1,700 iterations), the pooled mean moves by **+0.03 points of
  % correlation energy** (61.03 → 61.06). The per-set change is −0.46 to +0.59 points, with about as many molecules
  going down as up (41 up, 38 down, 8 tied).
* **The network still beats the converged label** on 78/87 molecules (`rl4L`, +5.42 points) and 79/87 (`rl4n29f`,
  +6.44). On the untouched n29 test set `rl4n29f` leads by +6.5 points and wins 24/24 molecules.
* **The norb-19 margin is the fragile one.** Against the converged label, the networks gain +2.1 to +2.6 points and win
  7–10 of 12. A second converged label for the same molecules, the stored 500-iteration label warm-started to
  convergence, averages 68.28, 1.8 points above the cold-start converged label. Against the best of the three
  converged variants per molecule, the networks lead by only +0.1 to +0.6 points (still winning 6–10 of 12).
* **Labels are noisy at the 1–2 point level per molecule.** L-BFGS on this objective is chaotic. Six runs with
  identical inputs and settings (21 molecules, norb 15–17) differ by a per-molecule SD of 1.6 points and a range of up
  to 12 points after convergence; the set mean varies with an SD of 0.4 points. Even the per-molecule best of six
  converged runs, picked with knowledge of the energy, stays 5.1 points below `rl4L`, which beats it on 16/21.
* **The t2 residual cannot pick the better label.** Among labels of the same molecule, a lower fit residual goes with a
  better LUCJ energy in 52% of pairs (all label variants, 87 molecules) and in 39% of pairs among the six equivalent
  runs. Choosing the lowest-residual run gives 61.7 points, below the six-run mean of 62.0. The networks have a
  *higher* residual than the converged label on 85/87 molecules and the better energy on 78/87.

## Setup

* **Fits.** `pretrain/followups/task3_fit_labels.py` is `generate_compressed_targets.py` for one config
  (`square_reg0.005`), with the same canonical init and optimizer (scipy L-BFGS-B, default tolerances). It uses
  `maxiter` 5000 instead of 500 and adds a passive callback that snapshots the parameters at iterations 500, 1000, …
  "Converged" means scipy's success flag (`CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL`); all 87 fits ended that
  way. Every process was pinned to one core of scai2 with `taskset`.
* **Floating-point trajectories.** The L-BFGS trajectory depends on the XLA summation order (the CPU affinity). Where
  this run reproduces the stored label bit-for-bit at the stored label's iteration, the 500-iteration snapshot is the
  stored label (9/9 at norb 15–16, 17/22 at 17–18, 20/20 at n29 val). Where it does not (5 molecules at 17–18, all of
  norb 19 and the n29 test set, whose labels were made under different thread settings), the run's own 500-iteration
  snapshot (`lab500`) was also evaluated. At norb 19 and the n29 test set, the stored label was additionally continued
  to convergence from a warm start (`warm_conv`, fresh L-BFGS memory). At norb 19 the run's own 500-iteration label was
  also continued this way (`warmown_conv`).
* **Energies.** Exact GPU energies (complex64) at norb ≤ 19 through `rl_queue/exact_shared`: 300 energies. MPS χ256 with
  the frozen v1 engine (`--tn-impl v1`, as for every n29 number in the docs) through `rl_queue/task3_tn256`: 97
  energies, on 3 workers on scai2 GPU 4. All runs had 0 errors.
* **Controls.** Labels that had already converged before iteration 500 have identical parameters in both runs and were
  re-evaluated: exact energies reproduce to 0.0000 mHa (8 molecules); the two TN re-evaluations differ by 0.16 and
  0.10 mHa (+0.04 and −0.02 points). TN evaluation noise is therefore negligible here.
* **Candidates.** `label` (stored, 500 iterations), `labconv` (this run, converged), `pre4` (pretrained, 4 recycles),
  `rl4L` (norb 15–18 RL, step 25), `rl4n29` / `rl4n29f` (+ n29 TN RL, steps 15 / 30). Selection caveat: `rl4L` was
  selected on the norb 17–18 val set and `rl4n29` on the n29 val set. The untouched sets are norb 19 and the n29 test
  set.

## 1. 500 iterations vs converged vs the networks

Mean % of the CCSD correlation energy:

| set (energies) | label, 500 it. | label, converged | pretrained ×4 | norb 15–18 RL | + n29 RL step 15 | + n29 RL step 30 |
|:--|--:|--:|--:|--:|--:|--:|
| norb 15–16 val (9, exact) | 62.0 | 62.6 | 59.4 | 68.2 | 69.7 | 70.3 |
| norb 17–18 val (22, exact) | 61.7 | 62.0 | 51.1 | 69.1 | 68.8 | 68.6 |
| norb 19 test (12, exact) | 66.7 | 66.5 | 58.1 | 68.6 | 69.0 | 68.9 |
| n29 val (20, MPS χ256) | 57.2 | 57.4 | 53.8 | 63.7 | 65.7 | 65.5 |
| n29 test (24, MPS χ256) | 60.4 | 60.0 | 54.5 | 64.7 | 66.3 | 66.5 |
| all sets pooled (87, mixed) | 61.0 | 61.1 | 54.5 | 66.5 | 67.5 | 67.5 |

Paired differences in points (wins = molecules where the candidate is higher):

| set | converged − 500-it. label | pretrained ×4 − converged | norb 15–18 RL − converged | step 15 − converged | step 30 − converged |
|:--|--:|--:|--:|--:|--:|
| norb 15–16 val | +0.59 (3 up, 4 down) | −3.18 (6/9) | +5.63 (8/9) | +7.10 (8/9) | +7.65 (8/9) |
| norb 17–18 val | +0.29 (10 up, 6 down) | −10.91 (2/22) | +7.11 (19/22) | +6.77 (18/22) | +6.58 (17/22) |
| norb 19 test | −0.17 (6 up, 6 down) | −8.37 (1/12) | +2.07 (10/12) | +2.56 (7/12) | +2.38 (10/12) |
| n29 val | +0.18 (10 up, 10 down) | −3.59 (4/20) | +6.31 (19/20) | +8.32 (20/20) | +8.07 (20/20) |
| n29 test | −0.46 (12 up, 12 down) | −5.45 (5/24) | +4.72 (22/24) | +6.33 (24/24) | +6.54 (24/24) |
| all sets pooled | +0.03 (41 up, 38 down) | −6.57 (18/87) | +5.42 (78/87) | +6.46 (77/87) | +6.44 (79/87) |

At norb 19 and the n29 test set, "500-it." is the stored label, which lies on a different floating-point trajectory
from the new fits. On the new fits' own trajectory, convergence raises the energy there (+0.79 and +0.40 points, table
below). The negative entries in the first column come from the trajectory switch, not from convergence.

### Fits

| set | converged (new) | iterations, median (range) | CPU per fit (median, 1 core) | stored labels converged |
|:--|--:|--:|--:|--:|
| norb 15–16 val | 9/9 | 599 (445–878) | 6 s | 2/9 |
| norb 17–18 val | 22/22 | 527 (401–971) | 7 s | 8/22 |
| norb 19 test | 12/12 | 626 (458–1,000) | 14 s | 1/12 |
| n29 val | 20/20 | 760 (533–1,135) | 94 s | 0/20 |
| n29 test | 24/24 | 793 (579–1,700) | 100 s | 1/24 |

### 500 → converged on the same trajectory

| set | fits past 500 | Δ t2 residual, mean | Δ % corr, mean / median (range) | better / worse | Spearman(Δresid, Δ% corr) |
|:--|--:|--:|--:|--:|--:|
| norb 15–16 val | 7/9 | −0.0057 | +0.76 / −0.14 (−3.07 … +10.25) | 3 / 4 | +0.14 |
| norb 17–18 val | 14/22 | −0.0014 | +0.34 / +0.08 (−3.39 … +8.30) | 9 / 5 | +0.35 |
| norb 19 test | 11/12 | −0.0023 | +0.79 / +0.58 (−2.53 … +3.19) | 7 / 4 | +0.11 |
| n29 val | 20/20 | −0.0046 | +0.18 / +0.02 (−2.20 … +4.02) | 10 / 10 | +0.07 |
| n29 test | 24/24 | −0.0061 | +0.40 / +0.34 (−5.31 … +7.68) | 17 / 7 | −0.29 |

The extra iterations lower the residual by 0.001–0.006 on average, a fraction of a percent of its value (0.25–0.34).
The energy moves by several points on individual molecules, in either direction. On average the change is +0.2 to
+0.8 points, inside the trajectory noise of the next section. A tighter tolerance (12 molecules, 2–5× more
iterations) moves the energy by +0.15 points on average (mean |Δ| 0.28, max 1.02), less than the trajectory noise. The default
tolerances are therefore tight enough for energy comparisons.

## 2. Which converged label? (norb 19 and n29 test)

| candidate | norb 19 (12) | vs network `rl4L` / `rl4n29` / `rl4n29f` (wins) | n29 test (24) | vs `rl4L` / `rl4n29` / `rl4n29f` (wins) |
|:--|--:|:--|--:|:--|
| stored label, 500 it. | 66.66 | +1.89 / +2.38 / +2.21 (9, 10, 10) | 60.42 | +4.26 / +5.87 / +6.08 (22, 24, 24) |
| this run, 500 it. | 65.77 | | 59.55 | |
| this run, converged (`labconv`) | 66.49 | +2.07 / +2.56 / +2.38 (10, 7, 10) | 59.96 | +4.72 / +6.33 / +6.54 (22, 24, 24) |
| stored label → converged (`warm_conv`) | 68.28 | +0.27 / +0.77 / +0.59 (10, 6, 9) | 59.99 | +4.69 / +6.30 / +6.51 (22, 24, 24) |
| this run's 500 → converged (`warmown_conv`) | 66.83 | +1.73 / +2.22 / +2.04 (10, 6, 9) | – | |
| best converged variant per molecule | 68.45 | +0.10 / +0.60 / +0.42 (10, 6, 9) | 60.86 | +3.81 / +5.42 / +5.64 (22, 24, 24) |

The "vs network" columns are network − label in points (positive means the network is better), with the number of
molecules the network wins.

* At norb 19 the three converged labels of the same molecules average 66.49, 66.83 and 68.28. Which one counts as "the
  converged label" moves the network's margin by 1.8 points. The stored-label warm start happens to land high. With 12
  molecules and a per-molecule label SD of 1.5–2.5 points, set means of single label runs are uncertain by about
  ±0.5–0.7 points.
* Measured against the per-molecule best of three converged labels, the networks keep a small lead at norb 19 (+0.1 to
  +0.6 points), and `rl4L` still wins 10/12 molecules. Its mean lead is small because two molecules lose by 3.2 and
  7.2 points (`C2H3NO_rxn3840_P`, `C2H3NO_rxn3841_P`).
* At n29 the choice does not matter: every converged variant sits 4–6.5 points below the networks, on 22–24 of 24
  molecules.

The trajectory noise between the stored labels and this run, both at 500 iterations, is: norb 17–18 (5 molecules)
mean |Δ| 0.76, max 2.69 points; norb 19 (12) mean |Δ| 1.32, max 4.38; n29 test (24) mean |Δ| 2.15, max 7.51. The stored
label and the cold-start converged label differ by mean |Δ| 1.18, 2.54 and 3.01 points on those sets.

## 3. Run-to-run spread of the label (norb 15–17, 21 molecules × 6 runs)

Identical inputs and settings (canonical init, `square_reg0.005`, default tolerances, maxiter 5000). The six runs are
the stored labels (multi-threaded XLA; a 2-thread run is bit-identical to them), a single-threaded run, and 4 runs from
the canonical init perturbed by relative 10⁻¹⁵ Gaussian noise. Exact energies, points of % correlation energy:

| quantity | 500 iterations | converged |
|:--|--:|--:|
| per-molecule SD over the 6 runs: mean / median / max | 1.54 / 1.02 / 6.27 | 1.62 / 1.66 / 4.76 |
| per-molecule range (best − worst run): mean / max | 3.78 / 16.09 | 4.02 / 11.79 |
| molecules where all 6 runs agree (range < 0.01) | 0/21 | 0/21 |
| set mean: stored run / min / max over runs | 61.48 / 60.72 / 62.29 | 61.60 / 61.53 / 62.38 |
| SD of the set mean over runs | 0.51 | 0.41 |
| set mean of the per-molecule best run (selected by energy) | 63.61 | 64.06 |

| network | mean % corr | vs mean converged label | P(network > a converged label) | > best of 6 converged | vs best of 6 (mean) |
|:--|--:|--:|--:|--:|--:|
| pretrained, 4 recycles | 55.14 | −6.87 | 0.25 | 3/21 | −8.92 |
| norb 15–18 RL (step 25) | 69.16 | +7.15 | 0.85 | 16/21 | +5.10 |
| + n29 TN RL (step 15) | 69.14 | +7.13 | 0.82 | 15/21 | +5.08 |
| + n29 TN RL (step 30) | 69.20 | +7.19 | 0.83 | 17/21 | +5.14 |

Convergence is not the issue: all these runs met the same convergence criterion. The objective has many local optima
with nearly the same residual but LUCJ energies several points apart, and floating-point noise decides which one a run
reaches. More iterations do not shrink the spread (per-molecule SD 1.54 → 1.62).

## 4. Residual vs energy

Computed by `pretrain/followups/task3_resid_energy.py` (from `summary.json` and the fit `.npz` files; output
`summary_resid_vs_energy.json` next to the energies).

| comparison (same molecule) | pairs / molecules | lower t2 residual ⇒ better energy | note |
|:--|--:|--:|:--|
| six equivalent runs, converged | 312 pairs, 21 molecules | 39% | picking the lowest-residual run: 61.73; six-run mean 62.02; highest-residual run 62.99; oracle 64.06 |
| six equivalent runs, 500 it. | 313 pairs, 21 molecules | 47% | lowest-residual 61.37; mean 61.44; highest-residual 62.50; oracle 63.61 |
| all label variants (stored, lab500, labconv, warm_conv, warmown_conv), 5 sets | 350 pairs, 87 molecules | 52% | Spearman(Δresid, Δ% corr) −0.10 |
| 500 → converged, same trajectory | 76 fits | – | Spearman −0.29 to +0.35 per set (table in section 1) |
| `rl4L` vs converged label | 87 molecules | – | network residual higher on 85/87 (mean +0.055), energy better on 78/87 (+5.42) |
| `rl4n29f` vs converged label | 78 molecules | – | residual higher on 76/78 (+0.057), energy better on 69/78 (+6.30) |
| `pre4` vs converged label | 87 molecules | – | residual higher on 69/87 (+0.020), energy better on 18/87 (−6.57) |

`rl4n29` / `rl4n29f` residuals are missing for the 9 norb 15–16 molecules, hence 78 there.

Per molecule, Spearman(resid, % corr) over the six converged runs averages +0.27 (median +0.37, negative on 8/21). A
lower residual tends to come with a *lower* energy. Across molecules, by contrast, the converged label's residual and
% corr are negatively correlated (Spearman −0.55, −0.66, −0.77, −0.31, −0.58 for the five sets). That is most likely a
molecule-difficulty effect (amplitudes that compress badly belong to molecules that 2-layer LUCJ correlates badly), not
a lever: within a molecule the residual carries no usable information about the LUCJ energy. This agrees with exp12 of
`docs/optimize_true_learnability_study.md`: the compressed-DF minima are not energy-equivalent, so the tie-breaker
should be energy, not residual.

## 5. What this changes

* The label baseline in `docs/rl_larger_molecules.md` and in HANDOFF §1 stands. Converged labels give the same set means
  to within ±0.6 points. The concern that the 500-iteration cap understates the labels is answered: it does not.
* The noise floor of a "label" baseline is about ±0.5 points on the set mean (n ≈ 20) and 1.5–2 points per molecule.
  Margins of +5 to +8 points (norb 15–18, n29) are far outside it. The +1.9 to +2.6 point norb-19 margin is not robust
  to the choice of label run: against the best of three converged runs it shrinks to +0.1 to +0.6. Task 7
  (`task7_rl_seeds.md`) finds a seed SD of 0.8 points for the same margin.
* A best-of-k label (k fits per molecule, selected by energy) is a stronger and fairer baseline than a single fit. Six
  fits raise the norb 15–17 set from 61.5–62.4 (one run) to 64.1, and the networks still lead by about 5 points.
  Selection by residual does not work (section 4), so a best-of-k baseline needs k energy evaluations per molecule.

## Files

* Code: `pretrain/followups/task3_fit_labels.py` (fits), `task3_run_fits.sh`, `task3_run_traj.sh`,
  `task3_run_noise.sh` (launchers), `task3_build_tasks.py` (energy-task pickles), `task3_analyze.py` (all tables;
  writes `summary.json`), `task3_report_tables.py` (the two headline tables), `task3_resid_energy.py` (section 4, added on resume),
  `task3_lists/` (name lists).
* Labels: `rhf_targets_compressed_conv/` (gitignored). `square_reg0.005/` holds the 87 converged fits with snapshots
  at 500, 1000, …; `tight/` the tighter-tolerance fits, `traj_{p1,p2,s1..s4}/` the six-run spread, `warm_n19/`,
  `warmown_n19/` and `warm_n29test/` the warm starts, and `own500_n19/` the norb-19 500-iteration snapshots as a
  labels dir.
* Energies: `pretrain/opt_true/results/followups/task3_converged_labels/energy_*.json`. The analysis output is
  `summary.json`; the residual-vs-energy statistics of section 4 are in `summary_resid_vs_energy.json`.
* Logs: `rl_runs/followups/task3_converged_labels/` (fit shards, queue evaluations).
* Reproduce the tables: `pretrain/.train_venv/bin/python3 pretrain/followups/task3_analyze.py` and
  `… task3_report_tables.py` (seconds, CPU).

## Resources used

* CPU: 273 label fits, 2.48 core-hours in total (the 87 main fits took 83 min), single-core processes on scai2.
* GPU: 300 exact energies on the shared exact queue (scai2 GPUs 6/7, scai7 GPU 1). 97 TN energies on 3 workers on scai2
  GPU 4 (00:43 → 04:02, ~3.3 h wall).
* All task-3 processes are stopped. The 3 TN workers on scai2 GPU 4 (`scai2-gpu4-a/b/c`) had been idle since 04:02 and
  were stopped at 09:47, which freed GPU 4. `rl_queue/task3_tn256` is empty.
