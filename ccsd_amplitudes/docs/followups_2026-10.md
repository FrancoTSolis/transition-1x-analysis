# Follow-up tasks of 5 Oct 2026: synthesis

*For Fang, 5 Oct 2026. Covers tasks 1–8 and the MPS sampler of HANDOFF §5. Detailed reports are in
`docs/followups/`, result files in `pretrain/opt_true/results/followups/<task>/`, logs in `rl_runs/followups/<task>/`.
The key numbers below were re-checked against those result files. Everything is measured unless it is marked
"extrapolated" or "estimated".*

Conventions: **%** is the share of the CCSD correlation energy that the LUCJ variational energy recovers, and
**points** are differences in that share. QSCI errors are in mHa. Tags follow HANDOFF §2. `label` is the stored
optimize=True label (500 L-BFGS iterations), `labconv` is the label run to convergence, and `cdf_lin` is Lin et al.'s
compressed-DF initialization (unregularized, her recipe for QSCI). Lin et al. is arXiv:2511.22476.

## 1. Bottom line

**What changed in our picture**

1. **The variational headline holds, but some margins are smaller than we quoted.**
   * Converged labels leave the baseline where it was: the pooled mean goes from 61.03 to 61.06 % over 87 molecules.
     The networks beat the converged labels on 78/87 molecules (`rl4L`, +5.4 points) and 79/87 (`rl4n29f`, +6.4).
     On the n29 test set the gain is +6.5 points, on 24/24 molecules. (Task 3)
   * At norb 17–18, the +7.4 we quoted reproduces across 3 seeds (+7.31 ± 0.15), but it is inflated by checkpoint
     selection. The honest size is **+5.7 ± 0.7** (late-training mean) or +5.8 ± 1.4 (step 50). (Task 7)
   * At norb 19 the margin is **about +2 ± 1**. It is +2.11 ± 0.82 at the val-selected step and +1.24 ± 1.60 at
     step 50. Against the per-molecule best of three converged labels it is only +0.1 to +0.6. (Tasks 3 and 7)
   * The policy is not at the optimum of its own objective. Per-molecule SPSA (500 evaluations) reaches 73.8 % from
     the label and 75.8 % from the policy, against 70.3 % for one network call. With 5,000 evaluations it reaches
     78.6–78.7 % on one molecule. NOMAD with Lin et al.'s settings reaches only 63.9 % from the label. (Task 2)
2. **A better variational energy does not give a better QSCI energy.**
   * Across the six parameter sets, the energy order and the QSCI order are unrelated (pooled Spearman −0.09).
   * Among the reasonable states, the QSCI error follows the number of distinct configurations sampled (Spearman
     −0.96 with strings per spin).
   * At a fixed subspace size, the label, pretrained and RL states are within 0.5–1.6 mHa of each other.
     (Tasks 1 and 2)
3. **A more accurate n29 reward (dm χ128) did not give a better policy.** This is confounded with a halved batch.
   χ64 already ranks GRPO groups correctly, and its gains carry over to χ256. (Task 6)
4. **The dm TN engine is accurate up to 44 orbitals**, and its cost grows more slowly with size than the cost model
   assumes. At 33–44 orbitals, χ128 is within 0.28–1.81 mHa of χ256. (Tasks 4 and 5)
5. **TN-sampled SQD works at 58 qubits.** For QSCI, χ512 MPS samples cannot be told apart from exact samples (bias
   0.00 ± 0.11 mHa at norb 15–16). The bottleneck is the diagonalization: at norb 29 each subspace holds
   1.0–2.0 M determinants. (Sampler)
6. **2-layer LUCJ barriers are poor.** Against CCSD(T), the barrier MAE is 14.5 kcal/mol for `rl4n29f` and 4.6 for
   CCSD. Every LUCJ candidate overestimates every barrier. (Task 8)

**Are we more competitive than Lin et al.?** On the variational energy, yes, by far. On QSCI, which is her metric, no.

| metric | ours | Lin et al.'s recipe | verdict |
|:--|:--|:--|:--|
| LUCJ energy, norb 15–16 (9 molecules) | one call (0.19 s): 70.3 % (`rl4n29f`) | `cdf_lin` −401 % (above HF). NOMAD, 500 evaluations from our label: 63.9 % (11–27 min per molecule) | far ahead |
| QSCI, 10⁵ samples, norb 15–16, vs FCI | `rl4L` 8.6 mHa (label 10.1, `rl4n29f` 10.0) | `cdf_lin` 8.0 mHa | level (error ratio 0.94×) |
| QSCI, 10⁶ samples, norb 15–16, vs FCI | `rl4L` 2.08 mHa | `cdf_lin` 0.69 mHa | behind (0.36×) |
| QSCI, 10⁵ samples, norb 17 (9 distinct), vs CCSD(T) | `rl4L` 9.2 mHa | `cdf_lin` 8.7 mHa | level |
| QSCI after per-molecule optimization of QSCI (her TN baseline with exact samples; 2 molecules, P / TS) | one call: 7.2 / 8.5 mHa | SPSA 3.2 / 4.1 mHa; NOMAD 6.0 / 7.9 mHa | behind by 4.0–4.4 mHa (SPSA) and 0.6–1.2 mHa (NOMAD) |

Our energy-trained RL improves the label's QSCI error by 1.2×. Her per-molecule TN optimization improves her start by
1.1–4.2× (median 2.4×). The comparison is loose: our ansatz is 2-layer square, hers 1-layer heavy-hex, on different
molecules.

**For the RL objective**

* **If the target is the variational energy**, keep the energy reward. There is 5.5–10.6 points of per-molecule
  headroom (task 2). At n29, keep the χ64 reward and spend the compute on batch size (task 6).
* **If the target is QSCI/SQD**, the energy reward is the wrong proxy, and the reward has to score QSCI. Options:
  * the QSCI energy at 10⁵ samples. That is one subspace per state, since the batches coincide. Measured parts:
    ~10 s of GPU plus 20–50 s on 4 CPU threads per reward at norb 15–17;
  * an energy reward plus a term for configuration coverage.

  Task 2 is a warning. The QSCI-optimal states win by spreading the samples, and their LUCJ energies collapse (P after
  NOMAD: −1037 %). At a fixed subspace size they are the worst of all parameter sets. A combined reward is therefore
  the natural first test. Beyond norb 17, a QSCI reward needs a fast SQD solver (SBD).

**For the NERSC budget (10,000 node-hours)**

* **No change to the request is needed.**
  * Task 4 measured 74 energies/h (TITAN Xp) and 126 (shared H200) at n29 χ128. That brackets the model's ~103 per
    A100.
  * Task 5 measured a throughput ratio of 2.36 from 29 to 44 orbitals; the model assumes 3.77. The RL line comes to
    2,390–2,920 node-hours at the model's anchor (2,800 budgeted) and 1,980–4,070 at the measured anchors.
  * The A100 anchor is still unmeasured and is the dominant uncertainty.
* **Possible saving.** χ64 rewards at 30–44 orbitals would save ~1,300 node-hours (budget note). Group ranking at χ64
  has not been tested at 33–44 orbitals, where χ64 sits 0.42–1.15 points from χ256.
* **Two cheap fixes decide between the low and the high end:**
  * let the 4 workers actually share the GPU (CUDA MPS daemon, or one process per GRPO group). Today they time-slice;
  * reuse the MPO.
* **SQD line.**
  * At norb 29 the measured subspaces have 1.0–2.0 M determinants, against ~4 × 10⁶ assumed.
  * χ512 sampling is validated.
  * The SBD cost is still unmeasured.
* **If the RL reward moves to QSCI**, the RL line has to be costed again. The model does not include that.

## 2. The tasks

### Task 1: noiseless SQD/QSCI with Lin et al.'s protocol

* **Question.** Does a better LUCJ energy give a better QSCI energy?
* **Setup.**
  * Exact GPU state vectors; 10⁵ samples in 10 batches × 4,000.
  * PySCF selected CI with a spin penalty (her `_sci` path). QSCI here is SQD without configuration recovery.
  * 6 candidates, 5 independent sample draws per state. Two sets: norb 15–16 (9 molecules, FCI reference) and
    norb 17 (9 distinct C2HNO molecules, CCSD(T) reference).
  * Extras: 10⁶ samples, and fixed-size subspaces (the top 100 / 300 strings per spin).
  * All 1,098 units are complete (2,043 diagonalizations). Expanse charged 378 core-hours.

| norb 15–16, vs FCI | LUCJ % | QSCI 10⁵ (mHa) | strings per spin | QSCI better than label | QSCI 10⁶ (mHa) |
|:--|--:|--:|--:|--:|--:|
| truncated CCSD | 16.6 | 85.0 | 83 | 0/9 | 48.2 |
| `cdf_lin` | −401.3 | 8.00 | 929 | 7/9 | 0.69 |
| label | 62.0 | 10.10 | 520 | – | 2.55 |
| `pre4` | 59.4 | 11.41 | 448 | 1/9 | 3.12 |
| `rl4L` | 68.2 | **8.61** | 611 | 8/9 | 2.08 |
| `rl4n29f` | 70.3 | 9.97 | 538 | 4/9 | 2.58 |

* **norb 17, vs CCSD(T), 10⁵ samples.**
  * `cdf_lin` 8.72 mHa, `rl4L` 9.23 (better than the label on 9/9), label 10.71, `rl4n29f` 12.14 (worse than the
    label on 9/9, although it is 18 mHa better variationally).
  * At 10⁶ samples: `cdf_lin` −1.61, `rl4L` 0.41, label 0.98, `rl4n29f` 1.65.
* **Caveats.**
  * 10⁵ samples contain only 1,400–3,900 distinct configurations. So all 10 batches are the same subspace (except for
    `cdf_lin`), and her min/max error bars collapse. The noise was measured from independent draws instead: SD
    0.27–0.48 mHa.
  * PySCF SCI was used instead of Dice (≤ 0.08 mHa effect, no change of order).
  * At norb 17, CCSD(T) lies at least 0.2–2.6 mHa above FCI, so the errors quoted there are too small.
  * One seed per policy.
* **Report.** [followups/task1_sqd_vs_lin.md](followups/task1_sqd_vs_lin.md)

### Task 2: per-molecule direct optimization (NOMAD, SPSA)

* **Question.** How far is one network call from the per-molecule optimum? What does Lin et al.'s per-molecule
  optimization buy?
* **Setup.**
  * NOMAD with her exact parameter lines (500 evaluations), with 1 OpenMP thread instead of 4.
  * SPSA, lightly tuned on one training molecule.
  * Objective A: the exact LUCJ energy, on 9 molecules at norb 15–16.
  * Objective B: QSCI of 10⁴ exact samples (her objective), on C2H3N_rxn2858 P and TS, starting from the label.
  * Final scores use the task-1 protocol (10⁵ samples).

| objective A (% , mean of 9 unless noted) | start | NOMAD 500 | SPSA 500 | SPSA 5,000 (1 molecule) |
|:--|--:|--:|--:|--:|
| from the label | 62.0 | 63.9 | 73.8 | 78.7 |
| from one call (`rl4n29f`) | 70.3 | +0.3 to +0.8 (3 molecules) | 75.8 | 78.6 |

| objective B, QSCI at 10⁵ samples (mHa vs FCI) | P | TS | LUCJ % (P / TS) | strings per spin |
|:--|--:|--:|--:|--:|
| label start | 8.0 | 9.0 | 61.8 / 56.7 | 449–509 (label and one call) |
| one call (`rl4n29f`) | 7.2 | 8.5 | 68.1 / 70.3 | |
| NOMAD on QSCI | 6.0 | 7.9 | −1037 / 21.7 | 559–1,123 (QSCI-optimized states) |
| SPSA on QSCI | 3.2 | 4.1 | −92 / −45 | |

* **Findings.**
  * One call is worth a median of 129 SPSA evaluations from the label, and more than 500 NOMAD evaluations on 8 of 9
    molecules.
  * The SPSA optima of objective A are 0.8–1.1 mHa worse in QSCI than their starts. They sample fewer distinct
    strings.
* **Cost.**
  * One call takes 0.19 s on one CPU core.
  * Objective A takes 6–28 min of a TITAN Xp per molecule; objective B takes 1.4–8.7 h.
  * That is 2 × 10³ to 2 × 10⁵ times one call.
* **Caveats.**
  * Each setting ran once (one seed).
  * NOMAD was untuned for 508–574 parameters.
  * Exact samples were used instead of χ = 50 MPS samples.
  * Two QSCI runs stopped early: NOMAD on P at 346/500 evaluations (killed, finalized from its checkpoint) and SPSA on
    TS at 478 evaluations (3 h cap).
  * Six NOMAD runs from the network start were not run.
* **Report.** [followups/task2_direct_opt.md](followups/task2_direct_opt.md)

### Task 3: optimize=True labels run to convergence

* **Question.** Does the 500-iteration cap understate the labels?
* **Setup.**
  * The same objective and initial guess, run with L-BFGS-B to scipy convergence. All 87 fits converged, after
    401–1,700 iterations.
  * Energies: exact at norb ≤ 19, v1 χ256 at n29.
  * Six identical runs on 21 molecules measured the run-to-run noise.

| set | label 500 it. | label converged | `rl4L` | `rl4n29f` | `rl4n29f` − converged (wins) |
|:--|--:|--:|--:|--:|--:|
| norb 15–16 val (9) | 62.0 | 62.6 | 68.2 | 70.3 | +7.65 (8/9) |
| norb 17–18 val (22) | 61.7 | 62.0 | 69.1 | 68.6 | +6.58 (17/22) |
| norb 19 test (12) | 66.7 | 66.5 | 68.6 | 68.9 | +2.38 (10/12) |
| n29 val (20, χ256) | 57.2 | 57.4 | 63.7 | 65.5 | +8.07 (20/20) |
| n29 test (24, χ256) | 60.4 | 60.0 | 64.7 | 66.5 | +6.54 (24/24) |
| pooled (87) | 61.0 | 61.1 | 66.5 | 67.5 | +6.44 (79/87) |

* **norb 19.** Three converged label variants average 66.49, 66.83 and 68.28. Against the per-molecule best of the
  three (68.45), the networks lead by only +0.10 to +0.60.
* **Label noise.** Identical runs differ by a per-molecule SD of 1.6 points (range up to 12), and the set mean has an
  SD of 0.4. The t2 residual does not predict which label is better: in 39–52 % of pairs the lower residual has the
  better energy.
* **Caveats.** `rl4L` was selected on norb 17–18 val and `rl4n29` on n29 val; only norb 19 and the n29 test set are
  untouched. The norb-19 set has 12 entries but only 9 distinct molecules (the norb-17 list likewise contains one
  C2HNO reactant 4 times).
* **Report.** [followups/task3_converged_labels.md](followups/task3_converged_labels.md)

### Task 4: TN reward throughput at n29 (the NERSC RL configuration)

* **Question.** How many dm χ128 energies per hour does one GPU deliver at 29 orbitals with 4 workers × 4 threads?
* **Setup.** 20 molecules × 2 candidates (label, `rl4L`) = 40 energies through the queue.

| GPU | wall (s) | energies per hour | GPU-s per energy |
|:--|--:|--:|--:|
| TITAN Xp, idle host (scai2) | 1,954 | 74 | 48.8 |
| H200, shared with another user's ~100 % job (scai7) | 1,147 | 126 | 28.7 |
| cost-model A100 (not measured) | – | ~103 | 34.8 |

* **Caveats.** One run per GPU. The A100 was not measured (scai3 was unreachable). Each task rebuilds its MPO.
* **Results and report.** `pretrain/opt_true/results/bench_dm_chi128_{titan,h200}_4x4.json`; the unit-cost table and
  caveats in [nersc_2027_gpu_budget.md](nersc_2027_gpu_budget.md).

### Task 5: the dm engine at 33–44 orbitals

* **Question.** Accuracy and cost at 33 / 37 / 44 orbitals and χ 64 / 128 / 256, and whether the cost model holds.
* **Setup.**
  * Three molecules (R, P, TS) per size, one formula per size, plus two 29-orbital anchors. Policy `rl4n29f`.
  * scai7 H200 (shared with another user's job): 66 runs, 0 failures.

| orbitals | χ64 − χ256 (mHa) | χ128 − χ256 (mHa) |
|--:|:--|:--|
| 29 | 2.47–2.49 | 0.62–0.69 |
| 33 | 1.92–5.19 | 0.39–1.44 |
| 37 | 1.95–5.02 | 0.28–1.81 |
| 44 | 2.39–4.42 | 0.61–1.35 |

| orbitals | energies per GPU-hour (4 × 4 queue, χ128, same shared H200) | ratio to 29 orbitals |
|--:|--:|:--|
| 29 | 126 (task 4) | 1.00 |
| 37 | 107 | 1.18 |
| 44 | 53 | 2.36 (model 3.77) |

* **Accuracy.** At 33–44 orbitals, χ128 is 0.06–0.41 points from χ256. At 44 orbitals (one TS), χ256 is 0.68 mHa
  above χ512. `discarded_sum` remains a usable error bar.
* **Scaling with size (one worker, χ128, 29 → 44 orbitals).**
  * MPS build ∝ n^2.07 (model 2.6);
  * ⟨H⟩ ∝ n^3.98 (model 4.0);
  * MPO build ∝ n^5.72 (68 s at 44 orbitals).
* **The 4 workers time-slice the GPU** instead of running concurrently.
* **RL line:** 2,390–2,920 node-hours at the model's anchor, 1,980–4,070 at the measured anchors, against 2,800
  budgeted.
* **Memory.** χ128 takes 2.4–6.1 GB per process and χ256 5.8–8.1 GB. Four χ256 workers fit a 40 GB A100.
* **Caveats.**
  * Measured on a shared H200, not an A100.
  * One formula per size; the 44-orbital C7H16 is an easy alkane.
  * The queue tests used 8 tasks per size.
  * The host CPUs (Zen 5) are faster than Perlmutter's.
* **Report.** [followups/task5_dm_scaling.md](followups/task5_dm_scaling.md)

### Task 6: n29 RL with dm χ128 rewards on one GPU

* **Question.** Does a more accurate reward give a better n29 policy than the Expanse χ64 run?
* **Setup.**
  * dm χ128 rewards from 4 workers × 4 threads on the shared scai7 H200.
  * 6 molecules × 8 per step for 24 steps: 144 GRPO groups, against 360 on Expanse. 8.0 h.
  * Evaluated with v1 χ256, as for every n29 number in the docs.

| % (v1 χ256) | label | start `rl4L` | Expanse χ64, step 15 / 30 | dm χ128, step 20 / 24 |
|:--|--:|--:|--:|--:|
| n29 val (20) | 57.2 | 63.7 | 65.7 / 65.5 | 64.7 / 62.9 |
| n29 test (24, untouched) | 60.4 | 64.7 | 66.3 / 66.5 | 65.7 / 63.7 |

* **Paired on test.**
  * Step 20 vs the start: +0.99 ± 0.30 (20/24 molecules better).
  * Step 20 vs Expanse step 15: −0.62 ± 0.21 (the Expanse policy is better on 17/24).
  * Step 24 vs the start: −0.98 ± 0.29.
* **Instability.** Val at χ128 went 63.98 → 56.84 (step 10) → 65.00 (step 20) → 63.19 (step 24).
* **Ranking agreement.** dm χ128 sits +0.28 to +0.35 points above v1 χ256 and ranks checkpoints the same way.
* **Cost.** Median 963 s per 48-energy step, i.e. 20 s per energy. The Expanse run took 262 s per 96-energy step on one
  node.
* **Caveats.** Batch size and reward are confounded, and there is one seed per run.
* **Report.** [followups/task6_n29_rl_dm.md](followups/task6_n29_rl_dm.md)

### Task 7: seed variation of the norb 15–18 exact-reward RL

* **Question.** Are the RL gains reproducible across seeds?
* **Setup.** The `grpo_large_rl4` recipe with seeds 1 and 2; seed 0 is the documented run. All energies are exact.

| gain over the label (points) | seed 0 | seed 1 | seed 2 | mean ± SD |
|:--|--:|--:|--:|--:|
| norb 17–18 val, val-selected step | +7.40 | +7.39 | +7.13 | +7.31 ± 0.15 |
| norb 17–18 val, late-training mean (steps 20–50) | – | – | – | +5.7 ± 0.7 |
| norb 17–18 val, step 50 | +5.88 | +4.27 | +7.13 | +5.76 ± 1.43 |
| norb 19 (untouched), val-selected step | +1.89 | +1.43 | +3.02 | +2.11 ± 0.82 |
| norb 19, step 50 | +0.78 | −0.07 | +3.02 | +1.24 ± 1.60 |
| norb 15–16 small val, val-selected step | +6.23 | +7.28 | +7.68 | +7.06 ± 0.75 |

* **Reading.** The val-selected +7.4 is the maximum of 11 noisy evaluations on the set it is reported on. The RL step
  itself is robust: every seed improves on its start policy at norb 19 (+0.9 to +4.0).
* **Caveat.** The small val set was the selection set of the start policy `rl4s`.
* **Report.** [followups/task7_rl_seeds.md](followups/task7_rl_seeds.md)

### Task 8: reaction energies and barriers (zero compute)

* **Question.** What do the existing LUCJ energies give for reaction energies and barriers?
* **Setup.**
  * The 10 complete R/P/TS reactions at norb 17–19, from existing exact energies.
  * Errors against CCSD(T), in kcal/mol, from `pretrain/rl/reaction_energies.py`.

| method | reaction-energy MAE | barrier MAE | max barrier error |
|:--|--:|--:|--:|
| Hartree–Fock | 14.2 | 20.8 | 70.3 |
| label | 15.4 | 17.0 | 38.2 |
| `rl4L` | 10.5 | 16.9 | 28.1 |
| `rl4n29f` | 10.0 | 14.5 | 23.9 |
| CCSD | 2.4 | 4.6 | 6.9 |

* **Findings.** For every LUCJ candidate the mean signed barrier error equals the MAE: every barrier is
  overestimated. 2-layer LUCJ recovers less correlation at TS geometries.
* **Caveat.** The 10 reactions are not independent. Four share one C2HNO reactant and four share one C2H3NO reactant.
* **Results.** `pretrain/opt_true/results/reaction_energies_n17_19.json`. The only write-up so far is HANDOFF §1.

### MPS sampler for TN-sampled SQD

* **Question.** Can we draw bitstrings from the dm-engine MPS for SQD at norb ≥ 20? In which basis are they? Do they
  match exact samples?
* **Built.**
  * `pretrain/rl/tn_sampler.py` does perfect sampling on the GPU. 10⁵ samples take 0.2–0.9 s at norb 15–16 and
    0.5–2.2 s at norb 29.
  * The unit tests are in `pretrain/rl/tests/test_tn_sampler.py` and pass.
* **Basis.** The bitstrings are occupations of the split-localized orbitals after the final orbital rotation. They
  are **not** MO occupations.

| QSCI shift vs exact sampling (norb 15–16, 4 molecules, `rl4n29f`, 7–14 sample sets each; mHa) | mean ± s.e. |
|:--|--:|
| MPS χ64 (subspaces 18–31 % smaller) | +1.65 ± 0.12 |
| MPS χ128 | +0.67 ± 0.13 |
| MPS χ256 | +0.20 ± 0.10 |
| MPS χ512 | +0.00 ± 0.11 |
| exact state sampled in the MO basis instead (−1.4 to +1.0 per molecule) | −0.36 ± 0.11 |

* **norb 29 (2 molecules).**
  * χ256 samples are within TVD 1.5–3.8 × 10⁻⁴ of χ512 samples.
  * 10⁵ samples hold 2,600–3,200 distinct configurations, which make a single 1.0–2.0 M-determinant subspace.
  * PySCF would need ~11 h on 2 cores per subspace (estimated), so no full subspace was diagonalized.
  * A truncated subspace (400 strings, 160 k determinants, 63 min) gives 86.7 % of the CCSD correlation energy,
    against 72.8 % for the LUCJ state. That is a demonstration, not Lin et al.'s setting.
* **Caveats.** 4 molecules and one policy. χ-basis QSCI is a valid variational estimator, but it is not her MO-basis
  protocol.
* **Report.** [followups/sampler.md](followups/sampler.md)

## 3. Next steps (ranked by value)

1. **Decide the target metric** (no compute): variational energy or QSCI. That choice fixes the RL reward (§1).
2. **Count distinct strings at n29.** For label, `cdf_lin` and `rl4n29f`, sample with `tn_sampler` at χ512 (≤ 20 s
   per 10⁶ samples). At norb 15–17 this count predicted the QSCI ranking (ρ −0.91 to −0.96).
3. **Prototype a QSCI-aware reward at norb 15–17,** in the small-set RL: energy plus a coverage term, or QSCI at 10⁵
   samples. From measured parts, one reward costs ~10 s of GPU plus 20–50 s on 4 CPU threads.
4. **Build and benchmark the GPU SQD solver (SBD).** It gates norb-29 QSCI (1–2 M determinants), a QSCI reward beyond
   norb 17, and the unit cost of the NERSC SQD line.
5. **Correct the claims in HANDOFF §1 and `docs/rl_larger_molecules.md`** (no compute):
   * norb 17–18: +5.7 ± 0.7 (late-training mean) or +5.8 ± 1.4 (step 50), not +7.4;
   * norb 19: about +2 ± 1;
   * converged labels leave the label baseline unchanged.
6. **NERSC measurements**, on Perlmutter or the scai3 A100:
   * the A100 anchor at 29 and 44 orbitals;
   * 4-worker throughput with and without the CUDA MPS daemon;
   * χ64 vs χ128 group ranking at 33–44 orbitals, which decides the ~1,300 node-hour saving;
   * molecule-affine task claiming for MPO reuse. This touches `reward_queue.py`, a core module, so it needs your OK.
7. **Separate batch size from reward precision at n29.** A v1 χ64 run with 6 molecules per step takes ~2 h on one
   Expanse node (extrapolated).
8. **Best-of-k labels at norb 19.** Each fit costs ≈ 14 s of CPU plus ~100 s of exact GPU energy. This settles the
   norb-19 margin.
9. **Task 2 leftovers:**
   * the 6 NOMAD runs from the network start (~2.5 h on scai2 GPU 3);
   * the QSCI objective from the network start;
   * repeat runs with other seeds.
10. **Rule for future RL claims:** run ≥ 3 seeds, and select checkpoints on a set other than the one reported, or
    report a fixed step or the late-training mean.

## 4. State of the tree and the machines (5 Oct, 13:25 PDT)

* **Nothing is committed.**
  * Untracked: `docs/followups/`, `pretrain/followups/`, `pretrain/opt_true/results/followups/`,
    `rl_runs/followups/`, `pretrain/rl/tn_sampler.py`, `pretrain/rl/tests/test_tn_sampler.py`, and this file.
  * Modified: `pretrain/rl/grpo_slot.py` (task 6's `--queue-kind` option; the default behaviour is unchanged) and
    `rl_runs/bench_dm_chi128_gpuutil.log`.
* **No follow-up process is running.**
  * scai2 GPUs 2–5 are idle.
  * scai2 GPUs 6/7 and scai7 GPU 1 run only the coordinator's `rl_queue/exact_shared` workers.
  * scai7 GPU 4 holds only another user's job.
  * On Expanse, the task-1 agent reports that no jobs are queued; the ControlMaster connection (pid 400977) closes on
    its timer.
