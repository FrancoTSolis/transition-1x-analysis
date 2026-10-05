# ML-LUCJ: GPU budget for the 2027 NERSC renewal

*Fang Sun, 4 Oct 2026. Unit: Perlmutter GPU node-hours (1 node = 4 × A100). The numbers come from a cost model built on our June–October 2026 timings (plain Python, available on request).*

## 1. Bottom line

**We recommend requesting 10,000 GPU node-hours (40,000 A100 GPU-hours) for ML-LUCJ in 2027.**
- The model's central estimate is 10,011 node-hours (80% range 7,266–12,734 over the unit-cost uncertainties).
- The two least certain unit costs are the GPU SQD solver and our tensor-network (TN) engine beyond 29 orbitals. We benchmark both first on Perlmutter, and would adjust test-set sizes if they come in higher.

**On the earlier 10,000.** Our measured unit costs support it. The updated breakdown keeps the original structure:
- Energy evaluation is the largest cost: RL rewards (~2,800) and the SQD/QSD tests in the style of Lin et al. (~3,000), which are now explicit lines.
- Training covers larger orbital-pair models with long schedules and ~160 sweep configurations (~2,480).

| Line in the original text | Original | Updated (rounded to 10) | Notes |
|---|--:|--:|---|
| Architecture development and pretraining | 2,500 | 1,200 | 10 production runs of ~300 epochs, models up to ~40M parameters |
| Architecture/hyperparameter sweeps | 2,000 | 1,280 | ~160 configurations × ~40 epochs |
| RL fine-tuning | 3,000 | 2,800 | Reward: LUCJ variational energy (exact ≤ 18 orbitals, TN beyond) |
| Inference sweeps, generalization benchmarks, timing | 1,000 | 3,230 | Now includes the Lin et al.-style SQD tests (1,330) and QSD tests (1,670 = 1.25 × SQD) |
| Restarts, failure recovery, ablations, production | 1,500 | 1,490 | 17.5% contingency |
| **Total** | **10,000** | **10,000** | |

Two refinements since the original text, both from our 2026 results:
- The RL reward is the LUCJ variational energy, computed exactly or with tensor networks; QSCI/SQD energies are used for evaluation. At affordable bond dimensions, SQD energies tracked sample diversity more than parameter quality.
- The inputs used so far are CCSD amplitudes; MP2/MRPT inputs remain an option.

## 2. Ready-to-paste "Planned GPU usage"

Lines are rounded to 10 node-hours and sum to the request.

> We request 10,000 Perlmutter GPU node-hours. Each Perlmutter GPU node contains four A100 GPUs, so this corresponds to 40,000 A100 GPU-hours. The estimate comes from a cost model built on our measured timings and scaled to the full Transition1x set (30,205 molecules, 15–44 orbitals, 30–88 qubits).
>
> Planned GPU usage:
>
> - 1,200 GPU node-hours for architecture development and production pretraining. The network trains on the compressed double-factorization (optimize=True) objective and teacher LUCJ parameters over 28,356 molecules. Our current 0.9M-parameter model takes about 0.2 A100-hours per epoch; production models (orbital-pair tokens, up to ~40M parameters) cost up to ~8× more per epoch. We plan 10 production runs (one-shot and recycled, square and heavy-hex outputs, input variants, seeds) of ~300 epochs.
> - 1,280 GPU node-hours for architecture and hyperparameter sweeps: about 160 configurations × ~40 full-data epochs, comparing orbital-token models, orbital-pair-token transformers, sparse attention variants, recycling, amplitude compression schemes and U/J output parameterizations.
> - 2,800 GPU node-hours for reinforcement-learning fine-tuning (GRPO) with the LUCJ variational energy as the reward. Rewards are exact GPU state-vector energies up to 18 orbitals and GPU tensor-network energies (MPS, bond dimension 128) from 19 to 44 orbitals. Training runs as a size curriculum (3 campaigns of 1,200 steps × 128 energies), plus a heavy-hex variant for hardware runs and short development runs: about 636,000 reward energies at 9–88 A100-seconds each.
> - 140 GPU node-hours for tensor-network energies (bond dimension 256) of four parameter sets (CCSD initialization, optimize=True compressed double factorization, network, network + RL) on all 1,849 validation and test molecules, plus model selection.
> - 1,330 GPU node-hours for sample-based quantum diagonalization (SQD) tests following Lin et al. (arXiv:2511.22476). These cover 1,200 states (300 held-out molecules × 4 parameter sets) sampled from GPU MPS simulations, with 10^5 samples, 10 batches of 4,000, and subspaces of about 4 × 10^6 determinants (cap 1.6 × 10^7) diagonalized on GPUs. They also cover configuration-recovery SQD on IBM quantum-processor samples for 27 circuits, and the noiseless emulation of those circuits.
> - 1,670 GPU node-hours for quantum subspace diagonalization (QSD) tests on the same sets, estimated at 1.25× the SQD cost (range 1.0–1.5×).
> - 90 GPU node-hours for inference sweeps, timing studies, a per-molecule optimization baseline, and porting and benchmarking on Perlmutter. The timing studies compare sub-second neural inference with per-molecule optimization, which takes 0.5–10 CPU-minutes for compressed double factorization and an estimated 2–100 GPU node-hours per molecule for tensor-network-optimized SQD.
> - 1,490 GPU node-hours (17.5% of the above) for checkpoint restarts, failure recovery, model ablations and final production runs.

**GPU readiness (for the separate ERCAP field).**
- **Runs on GPUs today** on our group's A100, H200, V100 and TITAN Xp nodes:
  - network training (PyTorch);
  - exact LUCJ energies (our complex64 GPU kernel: up to 18 orbitals on 40 GB A100s, 19 on 80 GB);
  - tensor-network energies (our PyTorch block-sparse MPS engine). The Hamiltonian expectation value runs in block2 on the node's CPU cores, with 4 workers × 4 threads per A100.
- **To be built in early 2027:**
  - a batched GPU sampler for the MPS (small new code);
  - the GPU SQD eigensolver SBD (open source; published GPU runs on Frontier and on A100s).
- **Fallback:** reduced SQD settings on CPU nodes.

## 3. Form: "Information needed"

**Input.**
- RHF-CCSD amplitudes for the STO-3G frozen-core valence space: 15–44 spatial orbitals (30–88 qubits; median 33 orbitals).
- Per molecule:
  - t1, shape (n_occ, n_virt): 56–484 values;
  - t2, shape (n_occ, n_occ, n_virt, n_virt): 3,136–234,256 values (median 70,756), i.e. 12 kB–0.94 MB in float32;
  - n_occ = 8–22 and n_virt = 7–22;
  - orbital energies, plus orbital frames derived from the geometry.
- Whole dataset: 2.3 × 10^9 amplitudes (9.2 GB) plus 0.24 GB of frames.
- Energy evaluation also reads the molecular integrals: ~300 GB for all molecules as stored today (dense float64), ~40 GB with 8-fold symmetry.

**Output.**
- A two-layer LUCJ circuit with square connectivity:
  - orbital rotations U, shape (2, n, n);
  - Jastrow matrices J, shape (2, n, n), restricted to (p, p+1) same-spin and (p, p) opposite-spin pairs.
- The final orbital rotation comes from t1 and is not predicted.
- In ffsim's parameterization: 2n² + 4n − 2 = 508–4,046 parameters per molecule (median 2,308).
- The network restricts U to real orbital rotations, so it emits n(n−1) rotation generators plus 2(2n−1) Jastrow values (926 at 29 orbitals).
- Model size: 0.92M parameters today; sweeps go up to ~40M.

**Overall.**
- 171 formulas (C/H/N/O, 3–7 heavy atoms).
- 10,073 reactions: all of Transition1x.
- 30,205 geometries (reactant, product, TS), of which about 21–22k are distinct.

## 4. What "full scale" means

**Data.**
- 30,205 molecules: the reactant, product and transition-state (TS) geometries of all 10,073 Transition1x reactions.
- 80% of molecules have 60 or more qubits. Only 96 (≤ 19 orbitals) fit exact state vectors.
- Formula-disjoint split: train 28,356 / val 675 / test 1,174. This is the official Transition1x split plus a 39–44-orbital holdout.

**Training.**
- Label-free pretraining on the compressed double-factorization (optimize=True) objective.
- GRPO RL with the LUCJ energy as reward, as a size curriculum in 3 campaigns: exact GPU energies up to 18 orbitals, TN energies at bond dimension χ128 for 19–44 orbitals.
- A heavy-hex variant of the model for IBM hardware.

**Testing** follows Wan-Hsuan Lin et al. (arXiv:2511.22476): 10⁵ samples per state in 10 batches of 4,000, and QSCI with ≤ 4,000 strings per spin (≤ 1.6 × 10⁷ determinants). On the QPU they use 10 runs × 10⁶ shots and 10 configuration-recovery iterations.

| Tier | Molecules × parameter sets | Method |
|---|---|---|
| TN energies | (1,174 test + 675 val) × 4 (CCSD init, optimize=True, NN, NN+RL) | ⟨H⟩ at χ256 |
| TN-sampled SQD | 300 test molecules (all 11 test formulas) × 4 = 1,200 states | MPS at χ512 (χ1024 from 38 orbitals), 10⁵ samples, d ≈ 4 × 10⁶ |
| QPU subset | 9 molecules (29 orbitals, 58 data qubits) × 3 = 27 circuits | 5 runs × 10⁶ shots; SQD with recovery; TN emulation |
| QSD | Same sets | 1.25× SQD (range 1.0–1.5×, my estimate); not measured |

**Where we deviate from Lin et al.**
1. **Sampling.** We sample from our MPS engine instead of exact state vectors.
   - Our near-half-filled valence spaces are beyond exact simulation from 20 orbitals: ≥ 10¹³ determinants at 50 qubits, versus 4 × 10⁹ for their 52-qubit N₂.
   - They used a χ = 50 MPS only inside per-molecule optimization.
2. **Basis.** Our MPS is compact only in split-localized orbitals, so samples and SQD are in that basis.
   - This is untested for SQD.
   - On the QPU, the basis change folds into the final orbital rotation.
3. **QPU runs.** We use 5 runs × 10⁶ shots per circuit, not 10. Configuration-recovery SQD runs on 3 of the 5 runs; pooling all runs into one SQD is the cheaper option.
4. **Circuits.** We plan 1-layer heavy-hex circuits on the QPU, as in Lin et al., whose circuits had two-qubit depth 110–156 at 52–60 qubits. All classical tests use our 2-layer square ansatz.
5. **References.** Beyond 19 orbitals we have no near-exact reference, unlike their FCI/SHCI.
   - We therefore report the fraction of the CCSD correlation energy and comparisons between parameter sets.
   - At these sizes the SQD/QSD tests mainly measure how SQD scales with system size and initialization quality.

**Included / not included.**
- **Included:** every line in §6.
- **Not included:**
  - DMC reference energies from the outline. Our RL minimizes the variational energy directly, so DMC would only be needed for reporting, and on GPUs it could be large.
  - NOQE downstream runs.
  - Larger basis sets.
  - RL with SQD energies as the reward (not planned; ~100× the RL cost).

## 5. Measured unit costs

| Workload | Hardware | Measured | Model value (A100) | Basis |
|---|---|---|---|---|
| Pretraining epoch, 30,034 molecules | Lab GPU (type not logged), shared hosts | 0.18–0.45 GPU-h (4 runs of the 0.92M model, 0.9–2.7 h each); 3.0M model with 2× batches: 0.08 | 0.20 A100-h | measured; A100 speedup 1.5× assumed (range 0.8–2.5×) |
| Exact LUCJ energy | A100 40 GB; H200 shared with another job | 3.0 s (17 orb.); 9–23 s (18); 43–57 s (19, H200) | 0.3–18 A100-s (15–18 orb.) | measured; 15–16 orb. derived (TITAN Xp ÷ 2.44) |
| TN energy, 29 orb., χ64 / 128 / 256 | 1 worker, TITAN Xp, host load 60–120 on 48 cores | 150 / 274 / 851 s; ⟨H⟩ 126 → 36 s on 4 threads; error 1.8 / 0.42 mHa at χ64 / 128 (15–18 orb.) | 35 / 108 A100-s at χ128 / 256 (4 workers per A100) | measured; A100 extrapolated |
| TN energy, older engine, χ256, 29 orb. | A100, 5 workers | 59–80 A100-s | check: model is 1.4–1.8× higher | measured |
| GRPO RL, 15–18 orb., 50 steps × 32 exact energies | 6 mixed GPUs | 1.5 h; ≤ 9.2 GPU-h | ≈ 4.8 A100-h | measured; A100 value derived |
| SQD diagonalization, PySCF, 29 orb. | 4 CPU threads | d ≈ 10⁵: 18–32 min; d ≈ 0.3–1 × 10⁶: 9–42 h (15 diag.) | 125 core-h at d = 10⁷ | measured; extrapolated ×100 in d |
| SQD diagonalization, GPU SBD (Walkup et al. 2026) | 8 × A100 | N₂, 3.1 × 10⁸ determinants: 309 s | 12 A100-min at d = 10⁷, 29 orb. (corners 3.4–32) | published; extrapolated |

Inference is negligible: 0.06–0.5 s per molecule on 4 CPU threads. Our 2026 pilot studies (≤ ~100 lab GPU-hours plus ~15–20 Expanse CPU node-hours) established these unit costs.

## 6. GPU budget by component (node-hours)

"Central" uses central inputs. P10–P90 come from 4,000 Monte-Carlo draws over all input ranges. They are per line and do not add up, and rows may not sum exactly because of rounding. Unit costs in the formulas are rounded to 0.1 A100-s; the Central column uses unrounded values, so a formula can differ from its row by 1.

| Component | Formula (central inputs) | Central | P10–P90 |
|---|---|--:|--:|
| Pretraining, production | 10 runs × 300 epochs × 8 (model size) × 0.20 A100-h ÷ 4 | 1,200 | 462–1,530 |
| Architecture/hyperparameter sweeps | 160 configs × 40 epochs × 4 × 0.20 A100-h ÷ 4 | 1,280 | 622–1,982 |
| RL, 15–19 orb. (exact to 18) | 3 × 150 steps × 16 molecules × 8 samples × 1.15 (eval) = 66,240 energies × 9.2 A100-s ÷ 0.75 ÷ 3600 ÷ 4 | 56 | 38–79 |
| RL, 20–29 orb. (TN χ128) | 3 × 300 steps: 132,480 energies × 28.6 A100-s ÷ 0.75 ÷ 3600 ÷ 4 | 351 | 190–664 |
| RL, 30–37 orb. | 3 × 450 steps: 198,720 × 54.6 A100-s ÷ 0.75 ÷ 3600 ÷ 4 | 1,004 | 540–1,813 |
| RL, 38–44 orb. | 3 × 300 steps: 132,480 × 87.7 A100-s ÷ 0.75 ÷ 3600 ÷ 4 | 1,075 | 577–1,949 |
| RL, heavy-hex variant | One 20–29-orbital stage: 44,160 × 28.6 A100-s ÷ 0.75 ÷ 3600 ÷ 4 | 117 | 60–245 |
| RL, development runs | 20 runs × 3,100 energies × 34.8 A100-s ÷ 0.75 ÷ 3600 ÷ 4 | 200 | 133–260 |
| TN energies (test/val/model selection) | 1,849 molecules × 4 at χ256 + 40 checkpoints × 100 at χ128 + 10% of test energies repeated at χ512; ÷ 0.85 | 140 | 117–173 |
| SQD, TN-sampled: MPS + samples | 1,200 states, MPS at χ512–1024, 10⁵ samples | 189 | 107–436 |
| SQD, TN-sampled: diagonalization | 1,200 × 10 batches at d = 4 × 10⁶ (GPU SBD); 0.1–4.2 A100-h per state | 478 | 130–877 |
| SQD of QPU samples + emulation | 27 circuits × 3 SQD runs × 10 iterations × 10 batches at d = 1.6 × 10⁷ (29 orb.) + TN emulation | 668 | 267–1,142 |
| QSD | 1.25 × the three SQD lines | 1,668 | 830–2,664 |
| Inference, timing, baselines | 30 checkpoints × 30,205 molecules; 300 timing runs; per-molecule TN optimization of 6 molecules (8.4 node-h each) | 54 | 30–76 |
| Porting and calibration | Build and benchmark SBD-GPU; TN engine on A100 | 40 | 27–58 |
| Contingency | 17.5% of the above | 1,491 | 998–1,891 |
| **Total** | | **10,011** | **7,266–12,734** |
| **Request (rounded)** | | **10,000** | |

**Monte-Carlo summary.** Median 9,575, mean 9,846, P10–P90 7,266–12,734 node-hours.

**TN cost per energy** (A100-seconds; n = number of orbitals):

cost(n) = [75.5·(n/29)^2.6 + 35.4·(n/29)^4 + 40·(n/29)^4 / 8] × 1.2 / 4

- **Terms**, measured at 29 orbitals and χ128:
  - 75.5 s: MPS build (151 s on a TITAN Xp, halved for an A100);
  - 35.4 s: ⟨H⟩ (block2 on 4 CPU threads);
  - 40 s: MPO, built once per 8-sample group.
- **Factors.** 1.2 is contention. The divisor 4 is workers per A100, because ⟨H⟩ needs 4 of the 16 CPU cores per GPU.
- **MPS exponent 2.6 (range 2.0–3.2).** This is the fit at χ128.
  - Fits: 2.0 (χ32, over 15–18 orbitals), then 2.2 (χ64), 2.6 (χ128) and 3.2 (χ256), each over 15–29 orbitals.
  - The χ256 fit is inflated, because the 15–17-orbital states are nearly exact at χ256 (discarded weight ~10⁻⁵).
  - Theory: n² at fixed χ.
- **⟨H⟩ exponent 4.0 (range 3.5–4.2).** The Hamiltonian MPO has bond dimension ∝ n² across n sites; we measured 3.6–4.2.
- **Result at χ128:** 11 / 35 / 52 / 75 / 131 A100-s at 20 / 29 / 33 / 37 / 44 orbitals.
- **Bond dimension.** At 29 orbitals the MPS build grows as χ^1.2 (χ64→128) and χ^1.8 (χ128→256). One χ512 run gives χ^1.3, probably too low; the asymptotic scaling is χ³. Above χ512 we assume χ^1.8 (range 1.3–3.0).

**SQD on GPUs** uses the SBD eigensolver.
- **Published port.** The GPU port of SBD (Walkup et al., arXiv:2601.16169; github.com/AMD-HPC/amd-sbd) reports a 63–95× speedup per Frontier node over CPU. Normalized to a Perlmutter GPU node (4 A100) against a 128-core CPU node, that is ~17–29×.
- **Our plug-in.** A separate Python wrapper (github.com/Qiskit/sbd-eigensolver-python) plugs into qiskit-addon-sqd as its `sci_solver`. Its A100 performance on our molecules is untested.
- **Per-determinant cost.** We take the geometric mean of the two published runs (N₂ on 8 A100s, H₂O on MI250X GPUs).
  - We scale it with the square of the number of single excitations per string (nlink), relative to N₂.
  - This scaling reconciles our measured CPU cost at 29 orbitals with the published GPU anchors and speedup: it implies 27× per node.
  - Linear scaling would give a total of 8,242, and nlink^2.5 would give 11,793.
  - The N₂ anchor's basis is not stated in the paper; we assume (10e,26o) cc-pVDZ.
- **Volume.** The central plan has 20,370 diagonalizations.

**Subspace size.** We assume 4 × 10⁶ determinants per batch (~2,000 strings per spin), with a range of 10⁶–1.6 × 10⁷.
- Lin et al.'s simulated subspaces stayed well below the cap: 242–929 strings per spin at the compressed initialization, and up to 2,136 after TN optimization.
- Our string spaces are far larger than theirs: 6.8 × 10⁷ strings per spin at 29 orbitals, versus 6.6 × 10⁴ for their N₂/cc-pVDZ. So we expect more unique strings, but not the cap.
- Noisy QPU samples fill the cap, so the QPU tier uses 1.6 × 10⁷.
- If every TN-sampled subspace also filled the cap, the total would be 13,843.

**Queues and node types.**
- Single-GPU jobs run in the shared QOS, which charges per GPU, so node-hours = A100-hours / 4.
- Two kinds of jobs go to the 80 GB nodes (same charge): χ256 energies at 37 or more orbitals, because four such workers may not fit in 40 GB; and the χ1024 MPS builds (76–88 qubits).

## 7. Outside the GPU request (rough)

**QPU (not a NERSC resource).**
- Shots: 27 circuits × 5 runs × 10⁶ shots = 1.35 × 10⁸.
- Time: ~5.9 QPU-minutes per 10⁶ shots (IBM's job-time estimate).
  - SQD needs ~800 minutes (P10–P90 510–1,250).
  - The same 1.25× factor adds ~1,000 minutes for QSD. Applying the QSD rule to QPU time is our extension, and Krylov-type QSD could need more.
  - Total ~1,800 minutes (~30 h).
- Price and access: ~$86k–172k at IBM list prices ($48–96 per minute), or ~4.5 default OLCF QCUP windows (400 minutes per 28 days).
- Connectivity: IBM Heron is heavy-hex. The heavy-hex model variant is in the RL line.

**CPU.**
- **SQD fallback.** Without the GPU solver, the same diagonalizations in PySCF would take ~31,000 CPU node-hours for SQD and ~39,000 for QSD, in place of ~2,600 GPU node-hours.
  - That is 27× per node, inside the published GPU speedup normalized to Perlmutter nodes (17–29×).
  - Capping max_dim at 2,000 on the QPU tier too gives ~18,000 + ~22,000. That is still not realistic, so the fallback would be smaller test sets and Dice.
  - About 800 CPU node-hours would cover SQD prototyping at reduced settings (1,000 states, max_dim 1,000).
- **RL rewards on CPU nodes.** Our older TN engine at χ64 runs ~1,300 energies per CPU node-hour, against ~410 per GPU node-hour for the χ128 GPU reward.
  - All TN RL rewards would then take ~640 CPU node-hours.
  - But the errors are 5–25× larger at 15–18 orbitals (2–10 mHa vs a 0.4 mHa median), and about 15 vs 0.5 mHa at 29 orbitals.
  - Matching the GPU reward's accuracy on CPU needs χ ≈ 256, roughly 6× the cost (rough estimate from our engine's measured χ scaling).
- **Labels and data.**
  - optimize=True labels for val + test: ~200 core-hours.
  - RHF/CCSD is done for all molecules (13.6 core-hours).
  - Hamiltonians take ~2 s each.
- **Reference energies.** Near-exact references beyond 19 orbitals (e.g. DMRG with block2 on the 300-molecule SQD subset) are CPU work not estimated here.
- **Already counted.** ⟨H⟩ for the GPU TN energies runs on the GPU nodes' own cores.

**Storage.**
- Amplitudes: 9.5 GB.
- Hamiltonians: ~300 GB dense (~40 GB with 8-fold symmetry).
- TN samples ~1 GB and QPU bitstrings ~2 GB.
- MPS states, if kept: up to ~0.7 TB.
- Checkpoints: tens of GB.

**Another route.** QIS@Perlmutter is a separate NERSC program for quantum-information projects. Its 2026 call offered up to 20,000 GPU node-hours per proposal, so it is worth checking for a 2027 call to cover the SQD/QSD part.

## 8. Assumptions and risks

- **GPU SQD solver (largest uncertainty).** We have not benchmarked it on our molecules.
  - How its cost grows with orbital count moves the total between ~8,200 and ~11,800. The subspace size moves it between ~9,100 and ~13,800.
  - Cost is linear in subspace size (max_dim²), batches and recovery iterations.
  - Our first Perlmutter task, budgeted under porting, is to benchmark SBD-GPU at 29–44 orbitals and d = 10⁶–1.6 × 10⁷, and to count unique strings in faithful MPS samples.
- **TN engine beyond our data.**
  - Timings exist only for 15–29 orbitals and χ ≤ 512 (a single χ512 run), from one TITAN Xp worker on a loaded host.
  - The A100 throughput is extrapolated, and checked against the older engine's A100 runs.
  - Four workers per A100 assume the MPS build is about half GPU work. If it is fully GPU-bound, rewards cost up to ~1.8× more.
  - Still untested: χ ≥ 1024 (memory), a batched GPU sampler (still to be written), and SQD in our localized-orbital basis.
  - A short A100 timing before submission (4 workers × 4 threads, 29 orbitals, χ128) would pin down the reward cost.
- **RL scope is a planning choice.** We plan 3 campaigns × 1,200 steps; our runs so far were 30–50 steps. A χ64 reward at 30–44 orbitals would lower the total by ~1,300; χ256 there would add ~4,800.
- **QSD.** The 1.0–1.5× rule is my estimate (Fang). We budget 1.25×; 1.0× gives 9,619 and 1.5× gives 10,403. The variant (NOQE-type or Krylov/SKQD) is still open; Krylov variants use deeper circuits, which need larger χ and more QPU time.
- **Lower or higher.** Fewer RL campaigns, pooled QPU SQD (one SQD per circuit saves ~1,200) or reduced SQD settings would lower the total; none is assumed. Larger basis sets (active spaces ≤ 44 orbitals; mostly CPU), DMC references or SQD-energy RL rewards would raise it.
