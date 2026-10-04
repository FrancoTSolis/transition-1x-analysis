# Energy RL on larger molecules: exact GPU rewards up to norb 19, tensor-network rewards at norb 29

*Follow-up to `chem_frame_pretraining.md` §6–7 (energy RL on norb 15–16), 3 Oct 2026.*

* **Code** (`pretrain/rl/`): `gpu_energy.py`, `tn_energy_v1.py` (the TN engine of every n29 number here),
  `tn_energy.py`, `reward_queue.py`, `launch_workers.sh`, `queue_eval.py`, `policy_dump.py`, `grpo_slot.py`,
  `report_tables.py`; tests in `pretrain/rl/tests/`; Expanse job scripts in `expanse/`.
* **Results**: `pretrain/opt_true/results/energy_{largeval,smallval,n19,n29val,n29test}_*.json`.
* **RL logs**: `rl_runs/grpo_large_rl4/`, `rl_runs/grpo_n29_tn/` (+ `.log`).
* All tables below are printed by `python3 -m pretrain.rl.report_tables`.

Percentages are of the CCSD correlation energy, $(E_\mathrm{HF}-E)/(E_\mathrm{HF}-E_\mathrm{CCSD})$. "Labels" are
optimize=True compressed double factorizations: L-BFGS, ≤ 500 iterations, λ = 0.005, the protocol of every label set
in this project.

## 0. Summary

| Question | Answer |
|:--|:--|
| Can the energy reward reach larger molecules? | **Exact, up to norb 19** (38 qubits), on one GPU. `gpu_energy.py` reproduces ffsim's LUCJ energy to ~1e-7 Ha in complex64. It takes 3–12 s at norb 17, 9–23 s at norb 18 and 54 s at norb 19 (45.7 GB state on an H200), 230–740× faster than one CPU core. Beyond that, an MPS at norb 29 (58 qubits): χ 64 for rewards, χ 256 for evaluation. |
| RL on norb 15–18 with exact rewards | Rewards came from 6 lab GPUs through a shared-file queue (1.4 h, 50 steps). On 22 held-out norb 17–18 molecules the score went from 65.2 to **69.1 %** (best step; 67.6 at the last), against **61.7** for the labels. No loss on the original norb 15–16 val. |
| Transfer to never-seen sizes | **norb 19**, exact (12 molecules never used for training or checkpoint choice): 68.6 % vs 66.7 for the labels (+1.9, better on 9/12). **norb 29**, MPS χ 256: val 63.7 vs 57.2 (+6.5, 19/20). On an untouched 24-molecule test set over 4 formulas: **64.7 vs 60.4** (+4.3, 22/24). The gap is the same at three truncation levels (+5.97 / +5.88 / +5.93), so it is not an MPS artefact. |
| RL at norb 29 with a TN reward | Expanse, 96 single-core χ-64 workers, 30 steps, 2.3 h. Val at χ 256 went from 63.7 to **65.7 %** (+8.5 over the labels). On the untouched test it went from 64.7 to **66.3–66.5 %**, better than the label on **24/24**, including the 2 formulas absent from RL. It also helped the smaller sizes: norb 15–16 68.2 → 69.7–70.3, norb 19 68.6 → 69.0, norb 17–18 unchanged within selection noise. |
| Throughout | One network call per recycle, 4 recycles (3 frozen pretrained + 1 trained). Each label is a separate L-BFGS fit per molecule. |

## 1. Exact LUCJ energies on one GPU (`gpu_energy.py`)

The reward of the earlier RL runs, an exact ffsim statevector on CPU, took 10–28 min per norb-17 energy on 16–25
CPU threads and could not run norb 18. `LUCJEnergyGPU` computes the same number,
`exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1))`, for closed-shell molecules on one GPU:

* **State.** The CI matrix $\psi$ (dim × dim, dim $=\binom{n}{k}$). For a spin-balanced UCJ state built on
  closed-shell Hartree–Fock it is exactly symmetric, and the engine relies on that throughout.
* **Orbital rotations.** The first one is applied to HF exactly, as a Slater determinant ($\psi = v v^T$ with $v$
  from $k\times k$ minors). Later rotations are factorized into Givens rotations (float64, on the CPU). Runs of
  rotations that touch only the lower orbitals go through one shared-memory numba kernel pass over contiguous
  column chunks. Then the engine transposes $\psi$ in place and applies the same pass again.
* **Diagonal Coulomb.** One phase per element.
* **Energy.** An eigen-factorized Hamiltonian, shifted by its exact HF expectation values, contracted tile by tile
  with batched GEMMs and a fused reduction. The HF shift brings complex64 from $10^{-5}$ to $10^{-7}$ Ha.

Validation (`pretrain/rl/tests/`):

| check | max \|ΔE\| (Ha) |
|:--|--:|
| 36 random systems (norb 4–10, indefinite ERIs, square / hex / all-to-all) + 6 STO-3G molecules vs ffsim | 7.5e-13 (complex128), 9e-7 (complex64) |
| independent verifier: 96 more random cases (norb 2–14, n_reps 1–3, every kernel path and fusion split) | 5.1e-12 (complex128); complex64 1.2e-6 on random indefinite Hamiltonians, ≤ 1.3e-7 on molecular integrals |
| all 390 stored norb 15–16 exact energies (every candidate type) | 1.0e-7 |
| 12 norb-17 (10,7) references | 1.0e-7 |
| full-size analytic Slater-determinant checks, norb 15–18 (incl. (9,9), dim 48620) | 1.6e-7 |
| same at norb 19 (11,8), dim 75582, 45.7 GB state (H200) | 4.6e-8 |

* The norb-19 run needed one fix: the phase kernel put the CI rows on CUDA grid axis y, which is capped at 65535.
* CPU float64 references (ffsim, pyscf) themselves carry ~1e-9 Ha of noise at norb ≥ 16: they move by 5e-10 under
  a global phase of the state. The engine is phase-invariant to the last bit.

Cost per energy (complex64; CPU ffsim on one core: 454 s at norb 15, 1660 s at norb 16, 1700–2700 s at norb 17):

| system (nocc, nvirt) | dim | TITAN Xp 12 GB | A100 40 GB | H200 (shared with another job) | peak memory |
|:--|--:|--:|--:|--:|--:|
| norb 15 (8,7) | 6435 | 0.7 s | | | 3 GB |
| norb 16 (9,7) / (8,8) | 11440 / 12870 | 2.2 / 2.9 s | | | 4.5 GB |
| norb 17 (10,7) / (9,8) | 19448 / 24310 | 7.3 / 11.5 s | 3.0 s / – | | 7 / 9 GB |
| norb 18 (11,7) | 31824 | 22.8 s | 9.2 s | | 12 GB |
| norb 18 (10,8) / (9,9) | 43758 / 48620 | – | 17.9 / 22.7 s | 8–10 / 12–14 s | 22 / 25 GB |
| norb 19 (11,8) | 75582 | – | – | 43–53 s | 55 GB |

## 2. Rewards from many machines: a shared-filesystem queue (`reward_queue.py`)

The GRPO driver (`grpo_slot.py --reward queue`) and the evaluator (`queue_eval.py`) write one pickle per energy
into `rl_queue/<name>/pending/` (file name `n<norb>_<na>_<nb>__<id>.pkl`). Workers on any machine claim tasks by
atomic rename into `claimed/<worker>/`, largest norb first and only up to the size their GPU can hold, and write
`done/<id>.pkl`. A task whose worker stops sending heartbeats is put back. Heartbeats come from a child process:
one TN energy can run 30 min, and block2 holds the GIL, so a heartbeat thread cannot run during it. TN collectors
also wait 30 min before they requeue a task. `launch_workers.sh <host> <gpu> <queue> <max_norb> ...` starts a worker
over ssh; `THREADS` sets its CPU threads. `WORKER_ARGS="--kind tn --chi 256 --tn-impl v1 [--zip-margin 3]"` starts a
tensor-network worker instead of an exact one.

**GPU limit.** On scai3–scai7 the worker refuses to start if user `fts` already has processes on 2 GPUs of that
machine. The count includes every job of the user, not only this queue. A worker re-checks every minute and exits
if the limit is exceeded. `python -m pretrain.rl.reward_queue guard` prints the current count.

## 3. RL on norb 15–18 with exact rewards

* **Starting point.** The 4-recycle policy of the small study (`rl_runs/grpo_slot_T4p`, step 30): 3 frozen
  pretrained recycles plus a trainable final recycle.
* **Train.** 51 molecules from 19 reactions (norb 15: 4, 16: 17, 17: 9, 18: 21; `pretrain/rl/large_train.txt`).
* **Val.** 22 molecules from 8 other reactions (norb 17: 12 C2HNO, norb 18: 10; `large_val.txt`). No val molecule
  is identical to a train molecule (checked by $E_\mathrm{HF}$).
* **GRPO.** As before: $G=8$, 4 molecules per step, $\sigma_K=0.01$, $\sigma_Z=0.005$ with mirrored pairs, one
  on-policy update, AdamW at lr 2e-5, 50 steps.
* **Rewards.** Exact GPU energies through the queue: scai3 A100 ×2, scai7 H200 ×2, scai2 TITAN Xp ×2 (norb ≤ 17).
  A step took 89 s on average (1.4 h for the run). On Expanse CPUs, one norb-17 energy had taken ~28 min.

| step | 0 | 5 | 10 | 15 | 20 | 25 | 30 | 35 | 40 | 45 | 50 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| val (22) mean % corr | 65.2 | 64.3 | 66.2 | 67.1 | 68.6 | **69.1** | 68.3 | 65.0 | 67.7 | 67.5 | 67.6 |
| val median | 64.0 | 66.7 | 66.3 | 66.7 | 68.6 | **69.3** | 68.2 | 65.1 | 67.8 | 67.4 | 67.4 |
| train (51) mean | 66.4 | | | | | | | | | | 69.0 |

Held-out norb 17–18 molecules (22, exact energies; "calls" = network evaluations per molecule):

| candidate | calls | mean % corr | median | min | vs label (mean) | better than label |
|:--|--:|--:|--:|--:|--:|--:|
| optimize=True label (L-BFGS ≤ 500 it.) | – | 61.7 | 61.7 | 46.0 | — | — |
| pretrained slot model, one shot | 1 | 45.9 | 49.7 | 9.7 | −15.8 | 1/22 |
| pretrained slot model, 4 recycles | 4 | 51.1 | 53.1 | 28.3 | −10.6 | 2/22 |
| small-set RL (norb 15–16), one shot | 1 | 63.4 | 64.1 | 49.1 | +1.7 | 13/22 |
| small-set RL (norb 15–16), 4 recycles | 4 | 65.2 | 64.0 | 59.1 | +3.5 | 15/22 |
| **norb 15–18 RL, step 25** | 4 | **69.1** | **69.3** | 61.2 | **+7.4** | **18/22** |
| norb 15–18 RL, step 50 (final) | 4 | 67.6 | 67.3 | 62.5 | +5.9 | 15/22 |

By size: the step-25 policy scores 69.9 % on the 12 norb-17 molecules and 68.2 % on the 10 norb-18 molecules. The
labels score 61.1 / 62.5 %.

* Step 25 is the best of 11 evaluations **on this same val set**, so 69.1 is optimistic. The final step (67.6)
  and the untouched test sets below are the fair numbers.
* The one-shot model trained the same way (`rl_runs/grpo_large_rl1`) got worse: val 63.4 → 64.0 → 62.5 → 60.3
  at steps 0/5/10/15, so it was stopped. The frozen recycles of the 4-recycle model give RL a better-conditioned
  start.
* **No forgetting at small sizes.** On the 9 val molecules of the original study (norb 15–16), the 4-recycle policy
  scores 67.4 % (median 66.0) before this run and 68.2 % (median 67.5) after it.

## 4. Size transfer with exact energies: norb 19

The only exactly solvable size above the RL sizes is norb 19 (38 qubits). The dataset holds 12 such molecules:
C2H3NO, 4 reactions, (11 occupied, 8 virtual). The 4 reactants are the same molecule, so there are 9 distinct ones.
None was used in any RL run or for picking a checkpoint.

* Hamiltonians and optimize=True labels were built for this test (`pretrain/rl/hamiltonian.py`,
  `generate_compressed_targets.py`, same settings as every other label set).
* Energies are exact, on one H200 that was shared with another user's job (55 GB, ~54 s per energy):

| candidate | calls | mean % corr | median | min | vs label (mean) | better than label |
|:--|--:|--:|--:|--:|--:|--:|
| optimize=True label | – | 66.7 | 67.7 | 54.6 | — | — |
| pretrained, one shot | 1 | 57.0 | 59.0 | 44.5 | −9.6 | 0/12 |
| pretrained, 4 recycles | 4 | 58.1 | 61.8 | 43.4 | −8.5 | 0/12 |
| small-set RL, one shot | 1 | 64.5 | 65.1 | 59.2 | −2.1 | 4/12 |
| small-set RL, 4 recycles | 4 | 65.7 | 67.3 | 58.5 | −1.0 | 7/12 |
| **norb 15–18 RL, step 25** | 4 | **68.6** | **69.9** | **61.8** | **+1.9** | **9/12** |
| norb 15–18 RL, step 50 | 4 | 67.4 | 68.7 | 61.9 | +0.8 | 9/12 |

* Counting the repeated reactant once (9 molecules), step 25 is +1.4 over the label (better on 6/9) and step 50 is
  +0.7.
* RL on norb 15–18 adds 2.9 points over the small-set RL policy at a size it never trained on.
* The margin over the labels is much smaller than at norb 17–18 or norb 29. The labels are strong on C2H3NO: 66.7 %
  here, against 61.7 % on the norb 17–18 val set and 57–60 % at norb 29. Two products are clearly worse than their
  labels (rxn3841_P 61.8 vs 69.1, rxn3840_P 68.5 vs 71.5).

## 5. norb 29 (58 qubits): tensor-network energies

### 5.1 The TN evaluator (`tn_energy_v1.py`)

`LUCJEnergyTN` gives the same LUCJ energy from an MPS that conserves $N_\alpha$ and $N_\beta$ separately. All numbers
in this section and below use the original zip-up engine, frozen as `tn_energy_v1.py`
(`reward_queue worker --kind tn --tn-impl v1`). `tn_energy.py` has since been replaced by a density-matrix engine
that is more accurate at equal χ (see its docstring and `pretrain/rl/tests/results/tn_v2/`); the evaluations here
were not re-run with it. As a control, three label energies recomputed with the frozen engine on other GPUs matched
the earlier values to ≤ 0.09 mHa.

* **Basis.** Split-localized orbitals: occupied and virtual MOs are localized separately (Boys from the geometry,
  or Edmiston–Ruedenberg from the integrals), and sites are put in Fiedler order. HF is an exact product state in
  this basis.
* **Jastrow factors.** Each square-mask term $e^{i z\, n'_A n'_B}$ couples two orbitals of the rotated frame. It is
  applied exactly as a bond-dimension-17 MPO (zip-up with an intermediate bond of margin × χ), then compressed to
  χ. The MPS is never orbital-rotated.
* **Final rotation.** The $t_1$ rotation and the basis change are folded into the Hamiltonian.
* **Energy.** A block2 expectation value; the operator is built once per molecule.

At norb 15, the exact final state needs a bond dimension of 159 per cut in the split basis (1e-5 discarded weight),
against 1123 in MO order. A Givens-circuit MPS in the network's chain order loses 29 % of the norm in a single
rotation. Without truncation it matches ffsim to 3e-13 Ha (complex128). Accuracy against exact energies, over 21
norb 15–17 tasks from the builder and an independent verifier:

| χ | typical error | worst error | Spearman vs exact inside 9-member GRPO groups |
|--:|--:|--:|--:|
| 32 | 5–20 mHa | 29.6 mHa (12.5 % corr) | 0.75–0.97 |
| 64 | 2–10 mHa | 16.6 mHa (7.0 %) | 0.93–0.98 |
| 128 | 0.6–4 mHa | 9.0 mHa (3.8 %) | 0.98–1.00 |
| 256 (zip margin 1.5) | 0.06–1.4 mHa | 2.6 mHa (1.1 %) | – |

* Errors are positive: truncation removes correlation. They grow with the correlation energy of the state.
* The zip margin is an accuracy knob like χ. On the hardest task at χ 64, margin 1.5 / 3 / 6 gives 16.6 / 9.5 /
  5.4 mHa.
* On n29 (C3H5N3), χ 32 / 64 / 128 / 256 recover 41.4 / 43.2 / 45.5 / 46.2 % of the correlation energy, and χ 256
  at margin 3 recovers 46.55 %. χ 256 is therefore not converged at n29: the residual is ≳0.5 % of the correlation
  energy.
* complex64 has a noise floor of ~0.1 mHa.
* Cost per n29 energy: χ 64 takes ~2.7 min on one CPU core (Expanse). χ 256 takes ~4 min per worker with 5 workers
  sharing one A100 or H200; 20 workers did 100 energies in 24 min. χ 256 at margin 3 is ~3× slower.

### 5.2 Size transfer to n29 at χ 256

Candidates are evaluated by TN workers through the queue (`queue_eval.py --kind tn`). The policies here never saw a
norb-29 molecule in RL; the pretrained slot model was trained on all sizes.

**Val set.** 20 molecules (9 C3H5N3, 11 C4H5NO; 20 reactions, none in the n29 RL train set). At χ 256 (zip margin
1.5):

| candidate | calls | mean % corr | median | min | vs label (mean) | better than label |
|:--|--:|--:|--:|--:|--:|--:|
| optimize=True label | – | 57.2 | 58.7 | 40.5 | — | — |
| pretrained, one shot | 1 | 50.5 | 50.4 | 36.3 | −6.7 | 1/20 |
| pretrained, 4 recycles | 4 | 53.8 | 52.6 | 47.3 | −3.4 | 5/20 |
| small-set RL (norb 15–16), 4 recycles | 4 | 61.2 | 60.1 | 54.3 | +4.0 | 15/20 |
| **norb 15–18 RL (step 25)** | 4 | **63.7** | **63.2** | **57.6** | **+6.5** | **19/20** |

**Truncation check.** χ 256 is not converged at n29, and truncation errors grow with correlation, so a candidate's
score could in principle depend on how MPS-friendly its state is. Label and RL were therefore re-evaluated at two
more accuracy levels:

| MPS | molecules | label | norb 15–18 RL | RL − label |
|:--|--:|--:|--:|--:|
| χ 128, margin 1.5 | 20 | 55.8 | 62.5 | +6.6 (19/20 better) |
| χ 256, margin 1.5 | 20 | 57.2 | 63.7 | +6.5 (19/20) |
| χ 128, margin 1.5 | 6 | 54.69 | 60.66 | +5.97 |
| χ 256, margin 1.5 | 6 | 56.19 | 62.07 | +5.88 |
| χ 256, margin 3 | 6 | 56.91 | 62.84 | +5.93 (6/6) |

* Going from margin 1.5 to 3 at χ 256 lowers the energy by 1.6–7.2 mHa per state (0.4–1.5 % of the correlation
  energy), about equally for label and RL.
* The gap moves by less than 0.1 point on every molecule. The ranking is a property of the states, not of the
  truncation. The absolute percentages are lower bounds, by roughly 1 % or more.
* Cost of margin 3 at χ 256: 11 min per energy on a V100 (with 2 CPU threads).

**Untouched n29 test set.** 24 molecules from 24 reactions not used anywhere, 6 per formula:
C4H5NO and C3H5N3 (16 occupied, 13 virtual) are the formulas of the n29 RL sets; C3H5NO2 (17, 12) and C5H5N
(15, 14) appear in no RL set. All 24 are distinct. Labels and Hamiltonians were built for this test
(`pretrain/rl/n29_test24.txt`). At χ 256:

| candidate | calls | mean % corr | median | min | vs label (mean) | better than label |
|:--|--:|--:|--:|--:|--:|--:|
| optimize=True label | – | 60.4 | 61.1 | 48.6 | — | — |
| pretrained, 4 recycles | 4 | 54.5 | 54.4 | 38.4 | −5.9 | 4/24 |
| small-set RL, 4 recycles | 4 | 62.5 | 62.6 | 54.8 | +2.1 | 19/24 |
| **norb 15–18 RL (step 25)** | 4 | **64.7** | **64.3** | 56.5 | **+4.3** | **22/24** |
| + n29 TN RL, step 15 (§5.3) | 4 | 66.3 | 66.0 | 58.7 | +5.9 | 24/24 |
| + n29 TN RL, step 30 (§5.3) | 4 | 66.5 | 66.7 | 59.5 | +6.1 | 24/24 |

| formula (nocc, nvirt) | label | pretrained ×4 | small-set RL | norb 15–18 RL | + n29 RL step 15 | + n29 RL step 30 |
|:--|--:|--:|--:|--:|--:|--:|
| C3H5N3 (16, 13) | 56.9 | 55.6 | 60.0 | 62.4 | 64.2 | 64.3 |
| C4H5NO (16, 13) | 58.4 | 51.9 | 61.3 | 63.8 | 65.9 | 65.7 |
| C3H5NO2 (17, 12), in no RL set | 62.5 | 47.7 | 62.4 | 65.6 | 66.4 | 66.6 |
| C5H5N (15, 14), in no RL set | 63.9 | 62.8 | 66.6 | 66.9 | 68.6 | 69.3 |

The policy RL-trained on norb 15–18 beats the per-molecule labels on all 4 formulas at norb 29. It never saw a
norb-29 molecule in RL, and two of the formulas appear in no RL set at all.

### 5.3 RL at n29 with a χ-64 reward (Expanse)

`expanse/grpo_n29_tn.slurm` runs on one 128-core Expanse node.

* **Rewards.** The policy driver (CPU, 8 threads) is served by 96 single-threaded TN reward workers: χ 64, zip
  margin 1.5, one cached evaluator each.
* **Start and data.** The norb 15–18 RL policy (step 25). Train: 40 n29 molecules (8 C3H5N3, 32 C4H5NO, 39
  reactions). Val: the 20 above.
* **GRPO.** $G=8$, 12 molecules per step (96 energies = one round of the workers), same exploration noise,
  lr 1e-5, 30 steps.
* **Cost.** 262 s per step, 2.3 h, no failed energy. Two earlier attempts died. In the first, the 112 workers
  inherited `OMP_NUM_THREADS=4` and oversubscribed the node. In the second, 3 cached evaluators per worker ran the
  node out of memory.

| step | 0 | 5 | 10 | 15 | 20 | 25 | 30 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| val (20) mean % corr, χ 64 | 60.75 | 60.77 | 60.47 | **62.78** | 62.04 | 61.31 | 62.55 |
| val median, χ 64 | 60.13 | 59.88 | 59.55 | 61.93 | 61.65 | 61.00 | 62.01 |
| train (40) mean, χ 64 | 58.5 | | | | | | 60.3 |

χ 64 underestimates these states by ~3 points (60.75 at χ 64 against 63.7 at χ 256 for the starting policy).
Re-evaluated at χ 256:

| MPS χ 256 | label | start: norb 15–18 RL | + n29 RL, step 15 (best χ-64 val) | + n29 RL, step 30 |
|:--|--:|--:|--:|--:|
| n29 val (20) | 57.2 | 63.7 | **65.7** (+8.5 vs label, 19/20) | 65.5 (+8.3, 20/20) |
| n29 test (24, untouched) | 60.4 | 64.7 | 66.3 (+5.9, 24/24) | **66.5** (+6.1, 24/24) |

* The χ-64 gain on val (+1.8 to +2.0) shows up unchanged at χ 256 (+1.8 to +2.0). The policy did not learn to
  exploit the χ-64 truncation.
* On the untouched test set the gain is +1.6 to +1.8, and it appears on every formula (table in §5.2). The step-30
  checkpoint was not picked by any criterion and is as good as step 15.

Exact energies at the smaller sizes, to check for forgetting:

| set (exact) | label | start: norb 15–18 RL | + n29 RL, step 15 | + n29 RL, step 30 |
|:--|--:|--:|--:|--:|
| norb 15–16, original small val (9) | 62.0 | 68.2 | 69.7 (9/9 above the start) | **70.3** (8/9) |
| norb 17–18 val (22) | 61.7 | 69.1 (selected on this set); 67.6 at step 50 | 68.8 | 68.6 |
| norb 19 test (12) | 66.7 | 68.6 | **69.0** (10/12 above the label) | 68.9 (10/12) |

Training at norb 29 cost nothing at the smaller sizes, and it helped norb 15–16 by 1.5–2.1 points. One policy, one
call per recycle, now beats the optimize=True labels from norb 15 to norb 29.

## 6. Caveats

* **n29 energies are MPS lower bounds.** χ 256 with zip margin 1.5 sits 0.4–1.5 points below convergence. The
  comparisons between candidates are stable to < 0.1 point across three truncation levels (§5.2).
* **Checkpoint choice.** Norb 15–18 step 25 was picked on the norb 17–18 val set, and n29 step 15 on the n29 val set
  (at χ 64). The norb-19 set and the n29 test set were never used for any choice. The step-50 and step-30 rows show
  the unpicked alternatives.
* **Labels.** The labels are the project's standard optimize=True protocol: L-BFGS capped at 500 iterations, and
  almost no fit reaches its tolerance. They minimize the $t_2$ compression residual, not the energy. "Beats the
  labels" is measured against this protocol.
* **Sets are small and some molecules repeat.** Each set has 9–24 molecules, and norb 17–19 is dominated by C2HNO and
  C2H3NO. The norb-17 val set holds the same C2HNO reactant 4 times, and the norb-19 set the same C2H3NO reactant 4
  times.
* **Run-to-run noise.** Each RL run used one seed. The norb 15–18 val curve wanders by ±2 points between evaluations,
  so differences of 1–2 points between single checkpoints are within that noise. The gains on the untouched sets
  (24/24 at n29) are not.
* **Numerics.** Exact energies are complex64 (~1e-7 Ha on molecules). The TN energies have a complex64 noise floor
  of ~0.1 mHa.
* **Shared GPUs.** The norb-19 energies ran on an H200 shared with another user's job (60 GB cap). Timings there
  include that contention.
* **Queue history.** Before the heartbeat moved into a child process, about 10 % of the long TN tasks were requeued
  and computed twice. That wasted compute; the collector keeps the first result for each task, so no result was
  affected. Six norb-18 tasks once ran out of GPU memory next to TN workers and were rerun on an H200.

## 7. Reproduce

```bash
# exact rewards / evaluation workers (respecting the 2-GPU limit on scai3-7)
bash pretrain/rl/launch_workers.sh scai3 1 $PWD/rl_queue/main 18
bash pretrain/rl/launch_workers.sh scai7 1 $PWD/rl_queue/main 19 complex64 60     # 60 GB cap on a shared H200
# RL on norb 15-18 (driver on CPU, rewards from the queue)
python3 -m pretrain.rl.grpo_slot --init runs_ot/slotall_T4/best.pt --prefix-T 3 \
    --init-policy rl_runs/grpo_slot_T4p/policy_best.pt --train-names pretrain/rl/large_train.txt \
    --val-names pretrain/rl/large_val.txt --reward queue --queue-root rl_queue/main --group-size 8 --batch-mols 4 \
    --sigma-k 0.01 --sigma-z 0.005 --antithetic --ppo-epochs 1 --lr 2e-5 --steps 50 --eval-every 5 \
    --out rl_runs/grpo_large_rl4
# candidates -> energies (exact or TN workers)
python3 -m pretrain.rl.policy_dump --names-file gauge_study/names_norb19.txt --policy init4:runs_ot/slotall_T4/best.pt \
    --tag pre4 --policy rl_runs/grpo_large_rl4/policy_best.pt --tag rl4L --labels-dir rhf_targets_compressed_n19 \
    --device cpu --out runs_ot/energy_tasks/n19_all.pkl
python3 -m pretrain.rl.queue_eval --dump runs_ot/energy_tasks/n19_all.pkl --queue-root rl_queue/main \
    --out pretrain/opt_true/results/energy_n19_all.json
THREADS=3 WORKER_ARGS="--kind tn --chi 256 --tn-impl v1" bash pretrain/rl/launch_workers.sh scai3 1 $PWD/rl_queue/tn256 29 complex64 "" 1 a
python3 -m pretrain.rl.queue_eval --kind tn --dump runs_ot/energy_tasks/n29val_base.pkl --queue-root rl_queue/tn256 \
    --out pretrain/opt_true/results/energy_n29val_base_chi256.json
# n29 RL on Expanse (one 128-core node, 96 single-threaded chi-64 TN reward workers; the v1 engine)
bash expanse/sync.sh
sbatch --export=ALL,RUN=grpo_n29_tn,INITPOL=rl_runs/grpo_large_rl4/policy_snap_for_n29.pt expanse/grpo_n29_tn.slurm
# its policies on every set (norb 15-16, 17-18, 19 exact; n29 val/test chi 256)
python3 -m pretrain.rl.policy_dump --names-file pretrain/rl/n29_test24.txt --policy rl_runs/grpo_n29_tn/policy_best.pt \
    --tag rl4n29 --policy rl_runs/grpo_n29_tn/policy_last.pt --tag rl4n29f --device cpu \
    --out runs_ot/energy_tasks/n29test_rl4n29.pkl
python3 -m pretrain.rl.report_tables          # every table of this document
```

norb-19 inputs: `python3 -m pretrain.rl.hamiltonian --names-file gauge_study/names_norb19.txt` (Hamiltonians);
`generate_compressed_targets.py --names-file gauge_study/names_norb19.txt --configs square_reg0.005 --out-dir
rhf_targets_compressed_n19` (labels; L-BFGS ≤ 500 iterations like every other label set).
