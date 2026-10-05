# MPS sampler for TN-sampled SQD at 29 orbitals (follow-up task "sampler", 5 Oct 2026)

**Question.** Lin et al.-style SQD/QSCI beyond exact state vectors (norb ≥ 20) needs bitstring samples from our
tensor-network state. Can we draw them from the dm-engine MPS of `pretrain/rl/tn_energy.py`, which orbital basis are
they in, and do they match exact samples once χ is large?

**Answer.**
* **Sampler.** `pretrain/rl/tn_sampler.py` draws perfect samples from the engine's MPS on the GPU: 10⁵ samples in
  0.2–0.9 s at norb 15–16 and 0.5–2.2 s at norb 29 (χ 64–512, one TITAN Xp). It also provides the integrals for SQD in
  the basis of the samples.
* **Basis.** The bitstrings are occupations of χ_j = Σ_p (F S)_pj φ_p, where φ are the MOs, F = exp(t1 − t1†) and S
  is the split-localized basis. They are **not** MO occupations (§2).
* **Distributions (norb 15–16).** On 4 molecules with policy `rl4n29f`, the χ 256 distribution is 1 − F =
  5·10⁻⁵–3·10⁻⁴ from the exact one, with total-variation distance 1–8·10⁻⁴. The shot-noise floor of 10⁵ samples
  is ~10⁻², so these samples are statistically indistinguishable from exact samples (§3).
* **QSCI energies (Lin et al.'s settings).** QSCI is more sensitive to the MPS truncation than the distribution is,
  because the subspace depends on which rare strings were drawn. Pooled over 7–14 independent 10⁵-sample sets per
  source and molecule, the MPS raises the QSCI energy above exact sampling by **1.65 ± 0.12 mHa at χ 64** (subspaces
  18–31 % smaller), 0.67 ± 0.13 at χ 128, **0.20 ± 0.10 at χ 256** (resolved only on the 2 of 4 molecules with the
  largest truncation, +0.3 and +0.5 mHa) and **0.00 ± 0.11 at χ 512**. One sample set scatters by 0.3–0.7 mHa
  (sd). Sampling the same exact state in the MO basis instead of the χ basis moves the QSCI energy by −1.4 to
  +1.0 mHa, either sign, with 1.3–2.5× larger subspaces (§4).
* **norb 29 (58 qubits).** The chain runs end to end on 2 molecules. χ 256 samples are within TVD 1.5–3.8·10⁻⁴ of
  χ 512 samples. 10⁵ samples hold only 2,700–3,200 distinct configurations, so all 10 of Lin et al.'s batches are the
  same 1.0–2.0 M-determinant subspace. PySCF's selected CI cannot diagonalize that within 2 h on ≤ 4 cores (estimated
  ≈ 11 h on 2 cores), so per the task I stopped at the subspace sizes. A truncated batch (the 400 most frequent
  strings, 160 k determinants) took 63 min and recovers 86.7 % of the CCSD correlation energy, against 72.8 % for the
  LUCJ state itself (§5).

## 1. What was built

| file | what |
|:--|:--|
| `pretrain/rl/tn_sampler.py` | `MPSSampler` (perfect sampling, amplitudes, probabilities); `build_mps` / `mps_energy` (the engine's state and its block2 energy, built once); `sampled_basis` / `sampled_basis_integrals` (the basis of the bitstrings, and H in it); bitstring helpers (`occ_to_ints` = ffsim `BitstringType.INT` / qiskit convention, `ffsim_sign`); exact references `LUCJStateGPU` (the gpu_energy engine with a basis change folded into the last orbital rotation) and `sample_ci_matrix`; SQD with Lin et al.'s settings (`LIN_SQD`, `sqd_batches`, `solve_batch`, `run_sqd`) |
| `pretrain/rl/tests/test_tn_sampler.py` | unit tests, all pass (`rl_runs/followups/sampler/test_tn_sampler.log`) |
| `pretrain/followups/sampler_validate.py` | norb 15–16: exact vs MPS samples (§3) |
| `pretrain/followups/sampler_seeds.py` | repeated independent sample sets for the QSCI statistics (§4) |
| `pretrain/followups/sampler_sqd.py` | QSCI stage (CPU) on saved sample sets |
| `pretrain/followups/sampler_n29.py`, `sampler_sci_bench.py`, `sampler_ccsdt.py` | norb 29: end-to-end run, selected-CI cost, CCSD(T) reference |
| `pretrain/followups/sampler_report.py` | prints the tables of this note |

Results are in `pretrain/opt_true/results/followups/sampler/` (`validate_*`, `sqd_*`, `n29_*`, `sci_bench_*`,
`ccsd_t_refs_sampler.json`, and `report.md`, the output of `sampler_report.py`). Samples (integer bitstrings, `samples/*.npz`) and logs are in
`rl_runs/followups/sampler/`. `tn_energy.py` and `gpu_energy.py` were **not modified**. The sampler uses only the
engine's public objects: `LUCJEnergySplitTN.state`, the MPS tensors `A`, the labels `q`, `move_center`, and `b2` for
the energy.

```python
from pretrain.rl.tn_energy import LUCJEnergyTN
from pretrain.rl import tn_sampler as S
ev = LUCJEnergyTN(h, g, const, norb, nelec, max_bond=256, device="cuda", name=name)   # dm engine, Boys/Fiedler basis
mps, Fr, info = S.build_mps(ev, U, Z, t1)              # Fr = final(t1) @ S: the sampled orbitals in MO coordinates
occ, logp = S.MPSSampler(mps).sample(100_000, seed=0)  # occ[i, j] = n_alpha + 2 n_beta of orbital j
ints = S.occ_to_ints(occ)                              # alpha = bits 0..n-1, beta = bits n..2n-1
hF, gF = S.sampled_basis_integrals(h, g, Fr)           # H in that basis (same constant)
batches = S.sqd_batches(ints, norb, nelec, seed=0)     # Lin et al.: 10 x 4,000, symmetrize_spin, max_dim 4,000
E = min(S.solve_batch(hF, gF, const, b, norb, nelec, spin_sq=0.0)[0] for b in batches)
```

**Algorithm.** Perfect sampling (Ferris & Vidal 2012). The constructor sweeps the block QR to the right end and back.
Sites 1…n−1 then become exact right isometries, even where the engine left zero-weight bond directions behind (the
CPU zip-up engine does this with cutoff 0). The conditional probability is then exactly p(s_i | s_<i) =
‖L_i A_i[s_i]‖² / ‖L_i‖². All samples of a batch advance one site at a time with dense GPU matmuls
(B, χ) @ (χ, 4χ) in complex64. The U(1)×U(1) labels keep every row of L_i inside one sector, so the dense product
wastes flops on zeros but needs no bookkeeping, and it is fast enough. Each draw carries its log-probability. N_α and
N_β are conserved exactly, so postselection never discards a sample. `amplitudes()` contracts given configurations;
the diagnostics below use it.

**Random numbers.** Uniforms come from numpy PCG64 on the host. torch's CUDA generator was not usable: with torch
2.4, consecutive float64 `torch.rand` calls on CUDA repeated parts of the stream, giving 24,576 duplicates in
4 × 65,536 × 6 draws. This biased the sampler's χ² test at 5σ; the same code on the CPU passed. Host uniforms also
make the draws identical on every device and for every batch size. Other GPU code in the project that relies on
float64 `torch.rand` on CUDA may want to check this.

## 2. Which orbital basis the bitstrings refer to

* **MO basis.** φ_p are the active RHF MOs of `rhf_hamiltonians/<name>.npz`. The LUCJ state of the energy jobs is
  |ψ⟩ = O(F U₂) e^{iJ₂} O(U₂†U₁) e^{iJ₁} O(U₁†)|HF⟩, where O(W) maps a†_i → Σ_j W_ji a†_j (ffsim) and F = exp(t1 − t1†)
  (`make_ucj_op(Z, U, "square", t1)`).
* **Basis of the MPS.** The engine never applies F or S to the MPS. Both are folded into its Hamiltonian,
  O(Fr)† H O(Fr) with Fr = F S. Here S holds the Boys-localized occupied and virtual MOs, localized separately, in
  Fiedler order; column j is site j. The MPS therefore holds c(x) = ⟨x; χ|ψ⟩ in the orbitals
  **χ_j = Σ_p (F S)_pj φ_p**, with AO coefficients `C_mo_active @ F @ S`. Equivalently,
  c = O(Fr)†|ψ⟩ = `ffsim.apply_orbital_rotation(psi, Fr.T, norb, nelec)`. Site j of a sample gives (n_jα, n_jβ) of χ_j.
* **Integrals for SQD.** SQD must use the same basis: h^χ = Frᵀ h Fr and (jk|lm)^χ = Σ Fr_pj Fr_qk Fr_rl Fr_sm (pq|rs),
  with the constant unchanged. They are real because t1 and S are real (`sampled_basis_integrals`).
* **Amplitude convention.** The MPS uses the interleaved JW order (0α, 0β, 1α, …) and ffsim the alpha-block order.
  For the same occupations, c_ffsim = phase · (−1)^{Σ_p n_β(p) Σ_{q>p} n_α(q)} · c_MPS. The probabilities are identical.
* **Verified** in the unit tests (n = 5–9, random Hamiltonians and split bases, CPU and CUDA engines):
  * MPS amplitudes match the ffsim χ-basis state to ≤ 2.5·10⁻¹³;
  * ⟨c|H^χ|c⟩ = E_MO to 4·10⁻¹⁴;
  * QSCI over all strings in the χ basis equals FCI in the MO basis to 7·10⁻¹³;
  * `LUCJStateGPU(basis=Fr)` reproduces c to 3·10⁻¹⁶ in complex128, and at norb 15–16 its complex64 χ-basis energy
    equals the MO energy to ≤ 1.3·10⁻⁶ Ha.
* **Not Lin et al.'s basis.** They measure the circuit's qubits, i.e. MO occupations, because their circuit includes
  the final orbital rotation. QSCI depends on the basis, so χ-basis QSCI is an equally variational but different
  estimator. MO-basis samples from this MPS would need O(Fr) applied to the MPS, a volume-law rotation ("What does not
  work" in the `tn_energy.py` docstring). §4 measures the MO-vs-χ difference with exact samples.

## 3. Validation at norb 15–16: MPS samples vs exact samples

**Setup.**
* **Molecules.** Four from `pretrain/rl/small_val.txt`:
  * C2H3N rxn2858 P and TS: norb 15, (8,8);
  * C2H4O rxn0724 P: norb 16, (9,9);
  * C3H4 rxn2391 TS: norb 16, (8,8).
* **Candidate.** `rl4n29f` (`runs_ot/energy_tasks/smallval_rl4n29.pkl`). Its LUCJ energies are 67–70 % of the CCSD
  correlation energy, 76–88 mHa above CCSD(T).
* **Exact samples.** χ-basis states from `LUCJStateGPU(basis=Fr)` in complex64, then inverse-CDF draws from |c|².
* **MPS samples.** `LUCJEnergyTN` dm engine, Boys/Fiedler basis from the shared cache `rl_runs/tn_basis_cache`,
  complex64, cutoff 1e-9, χ ∈ {64, 128, 256, 512}.
* 10⁵ samples per source; two exact sets, with seeds 1 and 2. One TITAN Xp (scai2 GPU 5), 4 CPU threads.

**Metrics.**
* **1 − F.** The fidelity |⟨ψ_exact|ψ_MPS⟩|², estimated without bias from the exact samples as
  ⟨ψ_m|ψ_e⟩ = E_{x∼p_e}[conj(ψ_m(x)/ψ_e(x))], with ± 2 standard errors. A second estimate from the MPS samples,
  E_{x∼p_m}[ψ_e/ψ_m], agrees within errors; it only works if the sign convention and the sampler are both right.
* **TVD full.** ½ Σ_x |p_e − p_m|, estimated as ½ E_{x∼p_e}|1 − p_m/p_e|. Standard errors are ≤ 3·10⁻⁴.
* **TVD top-10³.** The 1,000 most probable exact configurations plus one pooled "rest" bin.
  * "exact p" compares contracted MPS probabilities with exact ones, with no sampling noise.
  * "samples" compares the empirical MPS frequencies with p_e. The floor is the same comparison for a second exact
    sample set of the same size.
* **Distinct configurations, and distinct α∪β strings**, among the 10⁵ samples. The number of strings is √(QSCI
  subspace dimension).

| molecule (norb) | χ | E − E_exact (mHa) | discarded_sum | 1 − F | TVD full | TVD top-10³ (exact p) | TVD top-10³ (samples) | distinct configs | distinct α∪β | state build (s) | 10⁵ samples (s) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| C2H3N 2858_P (15) | exact | 0 | – | – | – | – | 0.0092 (2nd exact set) | 1359 / 1336 | 356 / 370 | – | – |
| | 64 | +1.386 | 8.2e-04 | 8.5e-04 ± 2e-04 | 1.7e-03 | 1.2e-03 | 0.0095 | 1313 | 294 | 7 | 0.19 |
| | 128 | +0.294 | 1.8e-04 | 2.9e-04 ± 8e-05 | 5.1e-04 | 2.6e-04 | 0.0100 | 1329 | 347 | 14 | 0.22 |
| | 256 | +0.047 | 3.7e-05 | 5.1e-05 ± 3e-05 | 1.6e-04 | 6.2e-05 | 0.0088 | 1389 | 361 | 45 | 0.37 |
| | 512 | +0.022 | 9.3e-06 | 4.0e-06 ± 1e-05 | 6.5e-05 | 2.4e-05 | 0.0086 | 1380 | 364 | 98 | 0.67 |
| C2H3N 2858_TS (15) | exact | 0 | – | – | – | – | 0.0108 (2nd exact set) | 1731 / 1796 | 438 / 439 | – | – |
| | 64 | +2.491 | 1.4e-03 | 1.5e-03 ± 2e-04 | 2.9e-03 | 2.2e-03 | 0.0116 | 1607 | 360 | 6 | 0.16 |
| | 128 | +0.741 | 3.7e-04 | 5.2e-04 ± 1e-04 | 8.9e-04 | 5.3e-04 | 0.0102 | 1733 | 423 | 14 | 0.21 |
| | 256 | +0.179 | 8.8e-05 | 1.9e-04 ± 6e-05 | 3.2e-04 | 1.2e-04 | 0.0101 | 1741 | 431 | 48 | 0.34 |
| | 512 | +0.052 | 1.7e-05 | 4.8e-05 ± 2e-05 | 1.1e-04 | 2.9e-05 | 0.0104 | 1759 | 447 | 116 | 0.94 |
| C2H4O 0724_P (16) | exact | 0 | – | – | – | – | 0.0074 (2nd exact set) | 836 / 896 | 250 / 279 | – | – |
| | 64 | +2.066 | 7.0e-04 | 2.7e-04 ± 1e-04 | 1.1e-03 | 1.1e-03 | 0.0071 | 836 | 243 | 7 | 0.17 |
| | 128 | +0.718 | 2.2e-04 | 8.9e-05 ± 7e-05 | 3.4e-04 | 3.0e-04 | 0.0073 | 869 | 260 | 15 | 0.22 |
| | 256 | +0.221 | 5.1e-05 | 5.2e-05 ± 4e-05 | 1.0e-04 | 7.2e-05 | 0.0065 | 855 | 262 | 44 | 0.39 |
| | 512 | +0.031 | 1.2e-05 | −7.0e-06 ± 9e-06 | 3.3e-05 | 2.4e-05 | 0.0067 | 855 | 272 | 114 | 0.72 |
| C3H4 2391_TS (16) | exact | 0 | – | – | – | – | 0.0129 (2nd exact set) | 2468 / 2463 | 571 / 587 | – | – |
| | 64 | +4.281 | 3.2e-03 | 3.2e-03 ± 4e-04 | 6.4e-03 | 4.6e-03 | 0.0153 | 2132 | 494 | 8 | 0.18 |
| | 128 | +1.218 | 9.6e-04 | 9.5e-04 ± 2e-04 | 2.2e-03 | 1.3e-03 | 0.0126 | 2374 | 556 | 16 | 0.23 |
| | 256 | +0.312 | 2.5e-04 | 3.4e-04 ± 1e-04 | 8.1e-04 | 3.3e-04 | 0.0133 | 2370 | 530 | 58 | 0.41 |
| | 512 | +0.082 | 4.9e-05 | 7.4e-05 ± 4e-05 | 2.1e-04 | 6.5e-05 | 0.0121 | 2449 | 570 | 201 | 0.79 |

* **The draws are exact draws of the MPS distribution.** The unit tests' χ² goodness-of-fit tests pass
  (p = 0.05–0.51). Each draw's log-probability equals log|amplitude|² of the contracted MPS to ≤ 3·10⁻⁵ (complex64;
  median 10⁻⁶). The canonical form's isometry error is ≤ 1.4·10⁻⁶.
* **1 − F tracks the engine's `discarded_sum`** within a factor ~2 and falls ~4× per doubling of χ. From χ 128 on,
  the TVD caused by the MPS (≤ 2·10⁻³) is far below the shot noise of 10⁵ samples (~10⁻²). The "samples" column
  sits at its floor.
* **The probability mass of the tail is where truncation acts.** Relative to exact, the MPS keeps 81–91 % of the mass
  in configurations of rank 10³–10⁴ at χ 64, 94–99 % at χ 128, 99–100 % at χ 256 and 100 % at χ 512; ranks 100–10³
  are ≥ 98 % already at χ 64. Correspondingly, χ 64 draws 0–14 % fewer distinct configurations and 3–18 % fewer
  distinct strings. From χ 256 on, the counts match the exact sets within their own scatter.

## 4. QSCI energies from MPS vs exact samples (Lin et al.'s settings)

**Settings.** `tn_sampler.LIN_SQD`, taken from her `scripts/quimb/*/lucj_compressed_t2.py`,
`src/lucj/quimb_task/lucj_sqd_quimb_task_sci.py` and `src/lucj/sqd_energy_task/lucj_compressed_t2_task_sci.py`:
* 10⁵ samples; postselection only (`max_iterations` 1, no configuration recovery);
* 10 batches of `samples_per_batch` = 4,000 distinct configurations, drawn without replacement by empirical frequency;
* α and β strings merged (`symmetrize_spin`) and truncated to `max_dim` = 4,000 by marginal count;
* PySCF selected CI in the product space (`solve_sci`, with `spin_sq` = 0.0 as in her state-vector task);
* the energy is the minimum over the batches.

`sqd_batches` builds the subspaces with the library's own first-iteration code (seed 0). A unit test checks that
min_b `solve_batch` equals `diagonalize_fermionic_hamiltonian` exactly. At norb 15–16 every sample set has fewer than
4,000 distinct configurations, so all 10 batches are the whole set. The QSCI energy is then a deterministic function
of which strings were drawn, and it scatters from one sample set to the next.

**One draw per source.** E_QSCI − E_CCSD(T) in mHa, with the subspace dimension in parentheses:

| molecule | E_LUCJ − E_CCSD(T) | exact χ #1 | exact χ #2 | MPS χ64 | MPS χ128 | MPS χ256 | MPS χ512 | exact MO |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| C2H3N 2858_P (15) | +75.90 | +7.305 (126,736) | +7.514 (136,900) | +8.851 (86,436) | +7.390 (120,409) | +7.270 (130,321) | +7.111 (132,496) | +6.224 (198,916) |
| C2H3N 2858_TS (15) | +87.87 | +6.741 (191,844) | +6.107 (192,721) | +9.226 (129,600) | +7.272 (178,929) | +7.208 (185,761) | +6.642 (199,809) | +4.996 (242,064) |
| C2H4O 0724_P (16) | +76.01 | +9.822 (62,500) | +7.260 (77,841) | +9.431 (59,049) | +9.077 (67,600) | +7.805 (68,644) | +8.202 (73,984) | +8.817 (195,364) |
| C3H4 2391_TS (16) | +84.02 | +6.040 (326,041) | +6.187 (344,569) | +8.423 (244,036) | +6.635 (309,136) | +6.934 (280,900) | +5.859 (324,900) | +6.835 (682,276) |

Two exact sets of the same state can differ by up to 2.6 mHa (C2H4O), so single draws cannot resolve sub-mHa
effects. **Pooled over independent sample sets** (all sets of `sqd_*{,_seeds,_seeds2,_seeds_mo}.json`, fresh seeds;
mean ± standard error of the mean, mHa vs CCSD(T), with the number of sets and the mean subspace dimension):

| molecule | exact χ | MPS χ64 | MPS χ128 | MPS χ256 | MPS χ512 | exact MO |
|:--|--:|--:|--:|--:|--:|--:|
| C2H3N 2858_P (15) | +7.45 ± 0.12 (n=14, 128k) | +8.39 ± 0.17 (n=7, 104k) | +7.50 ± 0.15 (n=7, 121k) | +7.22 ± 0.13 (n=13, 132k) | +7.20 ± 0.16 (n=7, 133k) | +6.08 ± 0.16 (n=7, 201k) |
| C2H3N 2858_TS (15) | +6.52 ± 0.14 (n=14, 197k) | +8.32 ± 0.29 (n=7, 150k) | +7.01 ± 0.17 (n=7, 181k) | +7.05 ± 0.12 (n=13, 186k) | +6.71 ± 0.11 (n=7, 193k) | +5.21 ± 0.23 (n=7, 247k) |
| C2H4O 0724_P (16) | +8.09 ± 0.20 (n=14, 77k) | +9.21 ± 0.20 (n=7, 56k) | +9.10 ± 0.27 (n=7, 66k) | +8.28 ± 0.14 (n=13, 73k) | +8.02 ± 0.17 (n=7, 76k) | +9.05 ± 0.18 (n=7, 195k) |
| C3H4 2391_TS (16) | +6.18 ± 0.09 (n=14, 324k) | +8.90 ± 0.11 (n=7, 224k) | +7.32 ± 0.20 (n=7, 276k) | +6.49 ± 0.10 (n=13, 297k) | +6.29 ± 0.16 (n=7, 314k) | +6.45 ± 0.08 (n=7, 673k) |

Difference from exact χ-basis sampling (mHa; mean ± s.e.; last row: mean over the molecules):

| molecule | MPS χ64 − exact χ | MPS χ128 − exact χ | MPS χ256 − exact χ | MPS χ512 − exact χ | exact MO − exact χ |
|:--|--:|--:|--:|--:|--:|
| C2H3N 2858_P (15) | +0.94 ± 0.21 | +0.05 ± 0.20 | -0.23 ± 0.18 | -0.25 ± 0.20 | -1.37 ± 0.20 |
| C2H3N 2858_TS (15) | +1.80 ± 0.32 | +0.49 ± 0.23 | +0.53 ± 0.19 | +0.19 ± 0.18 | -1.31 ± 0.27 |
| C2H4O 0724_P (16) | +1.12 ± 0.28 | +1.01 ± 0.34 | +0.18 ± 0.24 | -0.07 ± 0.26 | +0.96 ± 0.27 |
| C3H4 2391_TS (16) | +2.73 ± 0.15 | +1.14 ± 0.22 | +0.31 ± 0.14 | +0.12 ± 0.18 | +0.27 ± 0.12 |
| mean | +1.65 ± 0.12 | +0.67 ± 0.13 | +0.20 ± 0.10 | -0.00 ± 0.11 | -0.36 ± 0.11 |

Sources: exact χ = exact χ-basis state (`LUCJStateGPU(basis=Fr)`), sampled exactly; MPS χ = `tn_sampler` on the dm
engine's MPS; exact MO = the exact state in the MO basis (Lin et al.'s basis). Seeds are disjoint per source
(50,000 + r, 60,000 + 100 χ + r, 70,000 + r), so the sets are independent. Subspace dimensions are means over the
sets.

* **Truncation bias.** Averaged over the 4 molecules, the MPS raises the QSCI energy by 1.65 ± 0.12 mHa at χ 64,
  0.67 ± 0.13 at χ 128, 0.20 ± 0.10 at χ 256 and 0.00 ± 0.11 at χ 512. The bias falls 2.5–3.5× per doubling of
  χ. At χ 256 it is resolved on C2H3N TS (+0.53 ± 0.19) and C3H4 TS (+0.31 ± 0.14). These two have the largest
  `discarded_sum` at χ 256 (8.8·10⁻⁵ and 2.5·10⁻⁴, against 3.7–5.1·10⁻⁵ for the other two). At χ 512 no molecule
  differs from exact by more than 1.3 s.e.
* **Mechanism.** The bias follows the subspace size. Relative to exact sampling, χ 64 subspaces are 18–31 % smaller,
  χ 128 5–15 %, χ 256 −3 to +8 % and χ 512 −4 to +3 %. The MPS thins the tail of rare configurations (§3), and QSCI
  is built from exactly that tail.
* **Scatter.** One 10⁵-sample set of the exact state scatters by sd 0.34–0.74 mHa (range over 14 sets 1.3–2.6 mHa).
  A sub-mHa comparison between candidates or χ values therefore needs ≳ 5 sets per source, as here.
* **Basis.** With exact samples, the MO basis gives QSCI energies 1.4 and 1.3 mHa *lower* than the χ basis for the
  two C2H3N molecules, and 1.0 and 0.3 mHa *higher* for C2H4O and C3H4 (mean −0.36 ± 0.11). The MO subspaces are
  1.3–2.5× larger, because the χ-basis state is more concentrated. The χ basis is a different, equally variational
  QSCI estimator. It is not uniformly better or worse, and the difference (~1 mHa) is larger than the χ 256 truncation
  bias. TN-sampled QSCI at norb 29 can therefore be compared among candidates in the χ basis, but not with Lin et
  al.'s MO-basis numbers at better than ~1.5 mHa.
* **Practical rule (norb 15–16).** For QSCI within ~0.2 mHa of exact sampling, use χ 512, or χ 256 when
  `discarded_sum` ≲ 5·10⁻⁵. At norb 29 the χ 256 `discarded_sum` is 5.7·10⁻⁵ (C3H5N3) and 2.7·10⁻⁴ (C4H5NO). If
  the norb 15–16 relation carries over (not measured at norb 29), χ 256 QSCI there is biased by ~0.2–0.5 mHa.

The full output of `sampler_report.py` (all tables of this note) is saved as
`pretrain/opt_true/results/followups/sampler/report.md`.

## 5. norb 29 end to end (58 qubits)

**Setup.**
* Two molecules from the untouched `n29_test24.txt`: C3H5N3 rxn4472 P and C4H5NO rxn0158 P, both (16,16).
* Candidate `rl4n29f` (`runs_ot/energy_tasks/n29test_rl4n29.pkl`).
* No exact reference exists at this size. The χ 128 and χ 256 distributions are therefore checked against the χ 512
  MPS of the same state, with the estimators of §3 applied to the χ 512 samples.
* One TITAN Xp; 2 CPU threads for the host eigensolver and block2. Wall times are on a shared host (load ~30), and
  the C4H5NO χ 512 build overlapped with a CPU job.
* CCSD(T) for both molecules is from `sampler_ccsdt.py`, which uses the shared cache's code but writes its own file
  (`ccsd_t_refs_sampler.json`; about 1 s each on one thread).

| molecule | χ | E (Ha) | E − E_CCSD(T) (mHa) | % CCSD corr | discarded_sum | state (s) | 10⁵ / 10⁶ samples (s) | 1 − F vs χ512 | TVD vs χ512 | distinct configs | α∪β strings | QSCI subspace (each of 10 batches) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| C3H5N3 4472_P | 128 | -276.427739 | +141.1 | 72.59 | 2.9e-04 | 69 | 0.49 / 4.5 | 2.4e-04 ± 9e-05 | 6.9e-04 | 3,139 | 1,398 | 1,954,404 |
| C3H5N3 4472_P | 256 | -276.428490 | +140.4 | 72.76 | 5.7e-05 | 190 | 0.80 / 7.2 | 2.7e-05 ± 3e-05 | 1.5e-04 | 3,204 | 1,408 | 1,982,464 |
| C3H5N3 4472_P | 512 | -276.428612 | +140.2 | 72.78 | 2.1e-05 | 382 | 1.43 / 11.2 | (ref) | (ref) | 3,186 | 1,380 | 1,904,400 |
| C4H5NO 0158_P | 128 | -280.302109 | +134.4 | 67.55 | 8.9e-04 | 67 | 0.49 / 4.4 | 6.5e-04 ± 2e-04 | 1.4e-03 | 2,625 | 981 | 962,361 |
| C4H5NO 0158_P | 256 | -280.303517 | +133.0 | 67.91 | 2.7e-04 | 246 | 0.85 / 8.0 | 1.1e-04 ± 7e-05 | 3.8e-04 | 2,699 | 1,099 | 1,207,801 |
| C4H5NO 0158_P | 512 | -280.303934 | +132.6 | 68.02 | 9.5e-05 | 775 | 2.20 / 18.9 | (ref) | (ref) | 2,648 | 1,053 | 1,108,809 |

(1 − F and TVD against the χ 512 samples; the QSCI subspace is the same in all 10 batches.)

* **Energies.** The dm engine's χ 256 energy for C3H5N3 is 1.9 mHa below the v1 zip-up value in
  `energy_n29test_rl4n29_chi256.json` (−276.42654), consistent with §5.1 of `rl_larger_molecules.md`. Going from
  χ 256 to 512 lowers it by 0.12 mHa (C3H5N3) and 0.42 mHa (C4H5NO).
* **Sampling cost.** Sampling is cheap next to building the state: 10⁶ samples take 4.4–8.0 s at χ 128–256. The
  χ 256 samples are within TVD 1.5·10⁻⁴ (C3H5N3) and 3.8·10⁻⁴ (C4H5NO) of the χ 512 distribution, the same level as
  χ 256 against exact at norb 15–16.
* **Subspaces.** The states are concentrated: the HF_S configuration (the occupied split-localized orbitals, doubly
  occupied) has p = 0.81 (C3H5N3) and 0.83 (C4H5NO). 10⁵ samples hold 2,600–3,200 distinct configurations, fewer
  than 4,000, so all 10 of Lin et al.'s batches are the same set. That set has 981–1,408 α∪β strings, giving a
  0.96–2.0 M-determinant QSCI subspace.
* **QSCI cost.** One PySCF selected-CI σ build in these subspaces (`sampler_sci_bench.py`, 2 threads) takes 11 / 41 /
  158 s for 10 k / 40 k / 160 k determinants (the 100 / 200 / 400 most frequent strings), i.e. ∝ dim^0.96. For
  29 orbitals and 16 electrons per spin that is ~1 ms per determinant, because every (N−2)-electron intermediate
  string carries (n_virt+2)² = 225 links. The 160 k-determinant solve took 3,751 s, about 24 σ-equivalents (Davidson
  iterations plus RDMs). The full 1.98 M-determinant batch would need ~1,700 s per σ, about **11 h on 2 cores**
  (≳ 5–6 h on 4). That is well beyond the 2 h / 4-core budget, so I stopped there, as the task specified.
* **Truncated QSCI (not Lin et al.'s setting; demonstration only).** The 400 most frequent strings of the χ 256
  samples give E = −276.493491 Ha: 75.4 mHa above CCSD(T) and 86.7 % of the CCSD correlation energy, against
  140.4 mHa and 72.8 % for the LUCJ state. S² = 0.017; 160 k determinants; 63 min on 2 cores.

## 6. Caveats

* **Basis.** Every MPS-sampled QSCI energy is a χ-basis energy. These are variational and comparable among
  themselves, but not identical to Lin et al.'s MO-basis protocol. With exact samples at norb 15–16, MO-basis QSCI
  differs from χ-basis QSCI by −1.4 to +1.0 mHa per molecule (mean −0.4), and the MO subspaces are 1.3–2.5× larger
  (§4). The engine cannot do MO-basis TN-sampled SQD at norb 29.
* **QSCI is a tail statistic.** Its subspace consists of the strings seen at least once in 10⁵ samples, i.e. those
  with probability ~10⁻⁵ and up. An MPS can have a fidelity of 0.999 and still visibly thin that tail (χ 64).
  Increasing the shot count tightens the χ requirement further.
* **Scatter.** At 10⁵ shots a single sample set's QSCI energy scatters with sd 0.3–0.7 mHa, and two sets can differ
  by up to 2.6 mHa. Comparing candidates needs several sets per candidate, or more shots.
* **Scope.** One candidate (`rl4n29f`); 4 molecules at norb 15–16 and 2 at norb 29. More strongly correlated states
  (larger `discarded_sum`) need a larger χ for the same tail. `discarded_sum` is a usable a-priori flag, since it
  tracks 1 − F.
* complex64 is used throughout (states, sampler). Probabilities carry a relative error of ~10⁻⁶, far below every
  effect above.
* The norb-29 χ checks compare with χ 512, not with an exact state.
* `sqd_batches` calls a private function of qiskit-addon-sqd 0.13.1 (`_prepare_ci_strings`), so that its subspaces are
  identical to the library's. It refuses any other version.

## 7. What to run next

1. **A solver for the norb-29 subspaces.** Lin et al.'s batch is the full 1.0–2.0 M-determinant set here; PySCF's
   selected CI needs ~11 h on 2 cores for one (§5). Options: the GPU SQD solver (SBD; HANDOFF open task), or one
   Expanse node per subspace. That is an estimated ~1.5 h on 32 cores, ~50 core-hours, assuming ~50 % parallel
   efficiency of the σ build (not measured).
2. **Candidate comparison at norb 29 from subspace sizes first.** Task 1 found that, among reasonable LUCJ states,
   QSCI follows the number of distinct strings sampled (ρ ≈ −0.9 to −0.96 at norb 15–17). Sampling is cheap (≤ 20 s
   for 10⁶ samples at χ 512). Counting distinct strings per candidate (label, compressed DF, rl4n29f) at norb 29
   predicts the QSCI ranking before any diagonalization is paid for.
3. **χ choice.** Use χ 512, or χ 256 when `discarded_sum` ≲ 5·10⁻⁵, with ≥ 5 independent sample sets per
   candidate. At χ 512 the norb-29 state build takes 6–13 min on a TITAN Xp; sampling is negligible.
4. **MO-basis samples at norb 29** would need the final orbital rotation applied to the MPS, a volume-law
   operation. That is only worth building if an exact comparison with Lin et al.'s MO-basis protocol is required.
   Otherwise, χ-basis QSCI is a valid variational estimator in its own right.
