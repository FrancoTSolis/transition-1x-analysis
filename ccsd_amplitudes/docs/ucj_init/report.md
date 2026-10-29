# Network initialization for Avijit's UCJ optimizer

*29 Oct 2026, overnight study for Fang (to forward to Brandon and Avijit). Code: `ucj_init/`. Results:
`pretrain/opt_true/results/ucj_init/`. Logs: `rl_runs/ucj_init/`. Avijit's files (`avijit_simulator/`) are unmodified
(commit 92b76eb). "%" is the share of the CCSD correlation energy recovered (TS molecules) or of the FCI one (small
systems).*

**Task.** Make Avijit's UCJ simulator (`avijit_simulator/UCJ_optimizer.py`) converge better on our Transition1x
molecules. Compare a start from our network's LUCJ parameters with random and naive starts, and with Lin et al.'s
(Wan-Hsuan's) truncated (optimize=False) and compressed (optimize=True) starts.

## 0. Answer

1. **His simulator cannot converge from any start, ours included. The objective has no minimum.**
   * It minimizes the second-order BCH energy. For fixed orbitals that energy is *exactly* E_HF + gᵀM(U)g: a
     quadratic with no linear term.
   * Every stationary point therefore sits at E_HF.
   * Wherever correlation is possible, M has a negative eigenvalue, and the energy falls without bound.
   * This holds on all seven of his test systems (H4, LiH, C2), and on two TS molecules from five starts each.
   * A better start only runs away faster. From our network start, E_BCH2 − E_HF reaches −16.5 Eh after 25 cycles;
     the CCSD correlation energy is −0.26 Eh.
   * This is the "stop and report" case of `UCJ_context.md`: bug in §1, fix in §3.
2. **His ansatz cannot take our network output.**
   * It uses a real orbital rotation, one layer and a spin-free Jastrow. Projected onto it, every network or label
     start lands *above* Hartree–Fock (−28 to −55 %).
   * Exactly optimized, his ansatz holds 52–71 % of the FCI correlation on H4 (62 % on C2). A complex rotation in
     the same layer holds 90–99 %, and our 2-layer LUCJ 100 %.
   * The cause is the same as in point 1: with real rotations the first-order correlation is zero.
3. **The fix is the exact (variational) energy, in the LUCJ ansatz.**
   * Swapped into his own `RUCJ.kernel`, with nothing else changed, it reaches stationary points on H4 and LiH
     from every start. C2 needs more than his 100-cycle cap.
   * Two cheaper repairs fail:
     * a bounded polynomial completion of BCH2: its minima have the Hartree–Fock energy when evaluated exactly;
     * complex rotations: BCH2 stays indefinite at some rotations.
4. **On TS molecules, with the exact LUCJ energy and exact GPU gradients, our network is the best start we have.**
   * **Against Lin et al.'s starts it wins clearly** (16 molecules, norb 15–18):
     * it starts at 66–78 % before any step, and reaches 80 % in a median 10 gradient evaluations. Lin's truncated
       start needs 82–93, Lin's compressed start 82–98, random starts ~100;
     * at a fixed budget of 100 evaluations it is ahead of every Lin and random run;
     * it also ends higher: above the truncated start on 11/12 molecules (tied on 1), and above the compressed start
       on 8/8.
   * **With Avijit's own trust-region Newton driver** (§4.2) the gaps grow. Truncated and random starts need 87–103
     evaluations just to reach 70 %; the network start needs 4.
   * **Against our own optimize=True label the gain is real but small.**
     * One gradient step takes the label past the network's one-call value (2 evaluations).
     * At norb 15–16 the two trajectories then run together.
     * At norb 17–18 the network start stays ahead at a fixed budget (§4.1). Its final energy is higher on 9 of 12
       molecules (tied on 1, lower on 2).
   * **The zero start never moves**: it is a stationary point, the LUCJ version of his g = 0 problem.

**Where the network helps most:**
* on the larger molecules;
* with the trust-region Newton driver, where large early steps are expensive;
* whenever gradients are not available. Task 2 found one network call worth a median 129 SPSA evaluations.

## 1. Why the BCH2 objective cannot converge (the bug)

`_LocalModel.energy_elec` is

  E_BCH2(g, U) = ⟨0|H|0⟩ + ⟨x₀|H|x₀⟩ − ⟨h₀|xx₀⟩,  x₀ = X|0⟩, xx₀ = P X x₀, h₀ = H|0⟩,

with X = R̂(U) Ĝ R̂(U)† and Ĝ = ½ Σ g_pq n_p n_q. X is linear in g, so for fixed U:

* **No linear term.** For real H and real X, ⟨0|[H, X]|0⟩ = 0; the first-order BCH term vanishes identically.
* **Homogeneous quadratic.** E_BCH2 − E_HF = gᵀM(U)g exactly (checked to 1e-13).
* **No stationary point below E_HF.** ∂E/∂g = 2M(U)g = 0 forces E = E_HF.
* **Unbounded below whenever correlation is possible.** If M(U) has a negative eigenvalue, E_BCH2 → −∞ along it.
  A negative eigenvalue is exactly the condition for the UCJ state to lower the energy at small g.
* **It is the exact second-order Taylor polynomial of the exact UCJ energy.** It matches the exact (ffsim) energy of
  the same parameters at O(g²), to 6 digits. The exact energy is even in g; its bound comes from the O(g⁴) terms
  that BCH2 drops.

| system (all electrons) | E(g=0,U) − E_HF | dev. from t² scaling | negative modes of M at his seed | his default run, 50 cycles: E_BCH2 − E_HF | exact energy of that end point − E_HF | FCI corr. |
|:--|--:|--:|--:|--:|--:|--:|
| H4 linear 1.5 Å, STO-3G | 4e-16 | 3e-14 | 2 / 10 | −18.6 Eh | +0.41 Eh | −0.167 |
| H4 linear 2.0 Å | 4e-16 | 2e-14 | 2 / 10 | −60.7 | +0.13 | −0.322 |
| H4 square 1.5 Å | 0 | 1e-14 | 2 / 10 | −75.4 | +0.09 | −0.241 |
| LiH 1.6 Å | 0 | 1e-13 | 3 / 21 | −4.0 | +0.29 | −0.020 |
| LiH 3.0 Å | 2e-15 | 1e-13 | 3 / 21 | −23.6 | +0.07 | −0.088 |
| LiH 1.6 Å, 6-31G | 0 | 3e-14 | 11 / 66 | −0.82 | +0.32 | −0.019 |
| C2 1.25 Å | 7e-14 | 4e-14 | 4 / 55 | −0.44 | +0.18 | −0.221 |

None of these converged, and |g| grew every cycle to 10–30. In the exact energy every end point is *worse than
Hartree–Fock*.

On the TS-set molecules (frozen-core active space, norb 15, 25 cycles; CCSD correlation −0.264 Eh for TS and
−0.223 Eh for P):

| start | E_BCH2 − E_HF, TS | E_BCH2 − E_HF, P |
|:--|--:|--:|
| his default (CISD singles + relaxed-Hgg mode) | −1.4 Eh | −2.7 Eh |
| random | −2.1 | −0.9 |
| Lin truncated, projected (§2) | −13.2 | −1.6 |
| our label, projected | −11.1 | −3.1 |
| our network (rl4n29f), projected | −16.5 | 0.0 (†) |

(†) This one started above E_HF (+0.087 Eh). After 25 cycles it had reached the E_HF level, the only stationary
value BCH2 has, and was just tipping below it.

![BCH2 runaway](fig_bch2_runaway.png)

*Left: −(E_BCH2 − E_HF) per cycle from his default seed. The ticks at the right edge mark each system's CCSD
correlation energy. Right: the TS molecule from five starts; values above E_HF are drawn at 1e-6.*

**What does not fix it.**
* **Complex orbital rotations** (`bch2_complex_check.py`; one layer, random complex U).
  * They give a linear term, b ≠ 0, but M(U) stays indefinite at some rotations. Positive definite at 11 of 16
    random U (H4 ×2, LiH, C2); indefinite at the rest, where BCH2 is unbounded again.
  * Where it is bounded, the fixed-U minimum holds little: −0.003 to −0.010 Eh on C2 (FCI −0.221).
* **A bounded polynomial completion, "RQ2"** (`rq2_test.py`). This is the Rayleigh quotient of the truncated state
  (1 + iX − X²/2)|0⟩ in the CISD space:
  [E₀ + ⟨x₀|H|x₀⟩ − ⟨h₀|xx₀⟩ + ¼⟨xx₀|H|xx₀⟩] / [1 + ¼⟨xx₀|xx₀⟩].
  * It is variational in the CISD space, equals BCH2 to O(g²), and costs one more H application.
  * Its minima (47–99 % of the FCI correlation) sit at large g, where the truncated expansion means nothing.
  * The exact energy at those parameters is 0 ± 1 % on H4 ×3 and LiH ×2.
* **A norm cap on g.** The minimum then sits on the cap, whose radius is arbitrary.

## 2. The integration: our LUCJ output → his (u0, g0)

**Mapping, verified against ffsim.** His state is exp(iX)|HF⟩. In ffsim this is
`UCJOpSpinBalanced(n_reps=1, diag_coulomb_mats=[[g, g]], orbital_rotations=[R.T])`:
* ffsim's U is Rᵀ, because his `_Operator.rotate(R)` uses the rows of R;
* J_αα = J_αβ = g;
* the sign of the exponent does not matter for a real ansatz.

BCH2 equals the O(g²) coefficient of that ffsim energy on all seven systems (`bch2_reproduce.py`).

**Projection** (`bch2_ts.projected_seed`). Our output has two layers, complex U, spin-balanced J on the square
pattern, and a t1 final rotation. His ansatz has none of these. The best simple projection keeps:
* R = polar(Re U₁)ᵀ;
* g = Z₁ on the square pattern;
* no second layer and no final rotation.

It loses everything (`projection_loss.py`; exact energies, the 9 molecules at norb 15–16):

| start | as LUCJ | projected onto his ansatz | layer 1 only, complex U kept |
|:--|--:|--:|--:|
| network (rl4n29f) | 66–78 % | −28 to −44 % | 35–51 % |
| our label | 50–70 % | −28 to −55 % | 10–42 % |
| Lin truncated | 12–25 % | −1 to −2 % | 9–21 % |

LUCJ layers built from t₂ come in complex-conjugate pairs whose product gives real doubles amplitudes at first
order. A real rotation removes that term, as in §1.

**Capacity of each ansatz** (`ansatz_capacity.py`). Exact-energy L-BFGS from 4 random starts plus ffsim's truncated
start (2 + 1 for C2); best value, % of the FCI correlation:

| system | his ansatz: real U, 1 layer, spin-free | complex U, 1 layer | LUCJ, 2 layers, square (ours) | UCJ, 2 layers, dense J |
|:--|--:|--:|--:|--:|
| H4_lin_1.5 | 52.3 | 96.6 | 100.0 | 100.0 |
| H4_lin_2.0 | 71.0 | 99.1 | 100.0 | 100.0 |
| H4_sq_1.5 | 60.5 | 90.3 | 100.0 | 100.0 |
| LiH_1.6 | 91.6 | 99.6 | 100.0 | 100.0 |
| LiH_3.0 | 93.3 | 99.8 | 100.0 | 100.0 |
| C2_1.25 | 61.6 | 104.8 | (not run) | (not run) |

## 3. The fix, on his own driver: the exact energy

`ucj_init/exact_model.py` replaces only `_LocalModel` with the exact energy of the same state.
* It uses a dense state vector, autograd, and exact Hessian-vector products at κ = 0.
* `RUCJ.kernel` is unchanged: trust region, Steihaug-CG, accept/reject, his parameterization (accumulated U,
  occupied–virtual κ, packed g) and his convergence tests.
* Results are in `avijit_exact.py` → `avijit_exact.json`.

**Results**: % of FCI correlation reached / cycles used (cap 100). "(grad x)" marks runs stopped above
|grad|∞ = 1e-6. The starts:
* default: his own seed;
* zero: g = 0;
* random: two random seeds;
* truncated and compressed: ffsim optimize=False / optimize=True, one layer, projected as in §2.

| system | FCI corr. (Eh) | default | zero | random0 | random1 | truncated | compressed |
|:--|--:|--:|--:|--:|--:|--:|--:|
| H4_lin_1.5 | -0.167 | 19 % / 33 | -0 % / 0 | 37 % / 33 | 19 % / 40 | 38 % / 57 | 38 % / 58 |
| H4_lin_2.0 | -0.322 | 26 % / 15 | 0 % / 0 | 32 % / 25 | 22 % / 33 | 33 % / 34 | 69 % / 29 |
| H4_sq_1.5 | -0.241 | 33 % / 38 | -0 % / 0 | 38 % / 36 | 24 % / 28 | 23 % / 14 | 23 % / 11 |
| LiH_1.6 | -0.020 | 85 % / 68 | -0 % / 0 | 85 % / 100 | 85 % / 73 | 84 % / 32 | 91 % / 45 |
| LiH_3.0 | -0.088 | 63 % / 36 | -0 % / 0 | 72 % / 100 | 72 % / 100 | 72 % / 20 | 72 % / 100 |

* **H4 and LiH: every nonzero start ends at a stationary point**, with |grad|∞ 1e-8 to 1e-13 (one LiH run stopped at
  5e-7). His `converged` flag is False on some of them only because his triple test (step < 1e-6, |ΔE| < 1e-9 and
  gradient, all at once) stalls once the trust region cannot lower the energy any further.
* **C2 does not finish in 100 cycles** (|grad| 1e-3 to 1e-2, still descending).
* **The zero start "converges" at cycle 0 at E_HF.** g = 0 is a stationary point of the exact energy too.
* **The start decides the basin.** Example, H4 at 2.0 Å: his default start 26 %, random starts 22–32 %, the
  projected ffsim compressed start 69 %. Structured starts matter even in his ansatz.
* **Under BCH2 (same starts, same driver), 20 of 24 runs on H4 and LiH diverge**, down to −413 Eh. The other four
  are the zero starts, which "converge" at E_HF. In the exact energy, 19 of the 24 BCH2 end points are above E_HF.

## 4. TS molecules: exact LUCJ energy, exact gradients

**New engine:** `ucj_init/lucj_grad_gpu.py`, the exact energy *and gradient* of our LUCJ (adjoint method, GPU).
* The forward pass and inverse gates reuse the verified energy engine's Givens and phase kernels.
* H·ψ: the same-spin parts come from a dense string-space operator; the opposite-spin part from the energy engine's
  X tile with the full (P|Q) matrix, then one sparse product.
* The orbital gradient comes from spin-summed transition 1-RDMs, chained through U = expm(K) by the adjoint Fréchet
  derivative.
* Validation: against ffsim (complex128) the gradient is right to a relative 1e-6 (norb 8, 10, 12, central
  differences). At norb 15 the energy is within 1e-6 Ha of `gpu_energy`, and directional derivatives match its
  differences.

Cost per energy + gradient:

| norb | TITAN Xp | A100 | H200 |
|:--|--:|--:|--:|
| 15 | 2.1 s | | |
| 16 | 8.1 s | 3.9 s | 1.7 s |
| 17 | | 9.3 s | 4.1 s |
| 18 | | | 11.7 s |

**Optimizers** (`ucj_init/opt_exact.py`):
* L-BFGS-B, memory 20;
* Avijit's trust-region Newton logic, "trn":
  * radius 0.25, max 1.0; accept at ρ ≥ 0.1; ×0.25 / ×2 updates; Steihaug-CG, as in `RUCJ.kernel`;
  * Hessian-vector products by forward differences of exact gradients;
  * Steihaug capped at 20 inner steps (his cap is min(n_params, 50)).

Every gradient evaluation counts toward the budget, the Hessian-vector products included. The final orbital
rotation stays the t1 rotation for every start, as in task 2. Budgets:

| molecules | budget |
|:--|--:|
| TS, norb 15 | 400 |
| TS, norb 16 | 300 |
| TS, norb 17 | 250 (200 for rxn4442) |
| R/P, norb 16; norb 18 | 150 |

**Starts** (task-1 candidates: `cands_smallval.pkl`, `cands_n17c2hno.pkl`; norb 18: `ucj_init_n18.pkl`):

| tag | what | build time (norb 15) |
|:--|:--|--:|
| `zero` | U = I, Z = 0 | – |
| `random0`, `random1` | log U and J drawn from 0.1 N(0,1) | – |
| `truncated` | Lin et al.'s truncated t₂ start (ffsim `from_t_amplitudes`, optimize=False) | 2 s |
| `cdf_lin` | Lin et al.'s compressed DF (optimize=True, no regularization) | 42 s |
| `label` | our optimize=True label (regularized compressed DF) | 6 s |
| `pre4` | pretrained network (4 recycles) | 0.2 s |
| `rl4L`, `rl4n29f` | network + energy RL | 0.2 s |

### 4.1 L-BFGS

How to read the paired tables: each start is compared with `rl4n29f` on the molecules where both ran. "% at B" is
the best energy within B gradient evaluations; B = 0 is the start itself.

**All molecules** (norb 15–18):

| start | molecules | % at 0 | % at 10 | % at 25 | % at 50 | % at 100 | % at 150 | evals to 70 % (median) | evals to 80 % (median) | network ahead at 100 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| zero | 5 | -4.6 | – | – | – | – | – | never | never | 0/0 |
| random0 | 8 | -42.6 | 9.3 | 54.2 | 73.9 | 79.8 | 81.5 | 38 | 98 | 8/8 |
| random1 | 2 | -41.7 | 6.2 | 50.4 | 74.1 | 81.2 | 83.0 | 44 | 87 | 2/2 |
| truncated | 12 | 15.3 | 45.2 | 58.6 | 75.1 | 81.4 | – | 40 | 82 | 12/12 |
| cdf_lin | 8 | -350.6 | 62.0 | 72.3 | 77.0 | 79.7 | 81.1 | 20 | 98 | 8/8 |
| label | 14 | 59.2 | 80.1 | 82.2 | 83.3 | – | – | 2 | 8 | 9/12 |
| pre4 | 8 | 55.7 | 79.6 | 82.7 | 84.2 | 85.2 | 85.6 | 2 | 10 | 5/8 |
| rl4L | 5 | 66.7 | 79.8 | 82.5 | 83.8 | 85.0 | 85.4 | 2 | 11 | 2/5 |
| rl4n29f | 14 | 69.0 | 80.1 | 82.7 | 84.2 | – | – | 2 | 10 | – |

**TS molecules, norb 15–16:**

| start | molecules | % at 0 | % at 10 | % at 25 | % at 50 | % at 100 | % at 200 | evals to 70 % (median) | evals to 80 % (median) | network ahead at 100 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| zero | 4 | -5.3 | – | – | – | – | – | never | never | 0/0 |
| random0 | 4 | -41.2 | 12.6 | 55.0 | 74.0 | 79.6 | 81.9 | 34 | 129 | 4/4 |
| random1 | 1 | -38.7 | 6.7 | 54.9 | 78.2 | 84.0 | 85.0 | 38 | 56 | 1/1 |
| truncated | 4 | 16.1 | 43.2 | 57.8 | 74.6 | 80.8 | 83.9 | 42 | 105 | 4/4 |
| cdf_lin | 4 | -370.9 | 62.8 | 74.0 | 78.1 | 80.3 | 82.0 | 17 | 92 | 4/4 |
| label | 4 | 55.2 | 80.3 | 82.4 | 83.4 | 84.2 | 84.6 | 2 | 8 | 3/4 |
| pre4 | 4 | 58.3 | 80.0 | 82.8 | 84.2 | 85.4 | 85.9 | 2 | 10 | 2/4 |
| rl4L | 4 | 67.0 | 80.2 | 82.8 | 84.1 | 85.4 | 85.9 | 2 | 10 | 2/4 |
| rl4n29f | 4 | 68.6 | 80.2 | 82.9 | 84.2 | 85.4 | 85.9 | 2 | 10 | – |

**TS molecules, norb 17–18:**

| start | molecules | % at 0 | % at 10 | % at 25 | % at 50 | % at 100 | % at 150 | evals to 70 % (median) | evals to 80 % (median) | network ahead at 100 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| random0 | 3 | -46.6 | 6.7 | 53.3 | 73.2 | 79.8 | 81.8 | 42 | 105 | 3/3 |
| truncated | 3 | 11.5 | 46.1 | 57.9 | 72.3 | 79.5 | 81.9 | 43 | 84 | 3/3 |
| cdf_lin | 3 | -115.0 | 62.2 | 70.8 | 76.1 | 79.2 | 80.6 | 27 | 82 | 3/3 |
| label | 5 | 54.2 | 78.3 | 80.6 | 81.9 | – | – | 2 | 18 | 3/3 |
| pre4 | 3 | 49.8 | 79.5 | 83.2 | 84.8 | 85.7 | 85.9 | 3 | 9 | 2/3 |
| rl4n29f | 5 | 66.7 | 79.0 | 82.3 | 84.0 | – | – | 2 | 19 | – |

![exact optimization, L-BFGS](fig_exact_opt_lbfgs.png)

*Best-so-far % of the CCSD correlation energy against gradient evaluations (log scale beyond 10), TS molecules.
`pre4` and `rl4L` are omitted: after ~25 evaluations they lie on the `rl4n29f` curve.*

### 4.2 Avijit's trust-region Newton driver

Norb 15 (C2H3N P, TS) and norb 16 (C3H4 TS, C2H4O_rxn0724 TS), 400 evaluations:

| start | molecules | % at 0 | % at 10 | % at 25 | % at 50 | % at 100 | % at 200 | % at 300 | evals to 70 % (median) | evals to 80 % (median) | network ahead at 100 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| random0 | 4 | -38.3 | -2.3 | 12.1 | 37.6 | 67.6 | 77.8 | 79.8 | 103 | 281 | 4/4 |
| truncated | 4 | 15.7 | 39.8 | 45.5 | 57.8 | 72.2 | 80.0 | 82.1 | 87 | 207 | 4/4 |
| label | 4 | 58.3 | 78.7 | 81.7 | 82.6 | 83.8 | 84.5 | 84.7 | 3 | 18 | 1/4 |
| rl4n29f | 4 | 68.3 | 76.5 | 80.1 | 81.4 | 83.3 | 84.5 | 85.1 | 4 | 22 | – |

![exact optimization, trust region](fig_exact_opt_trn.png)

**A trust-region Newton method depends more on its start than L-BFGS does.** Its steps are capped at radius 1, and
each step costs up to 21 gradient evaluations. From Lin's truncated start it needs ~2× the evaluations L-BFGS
needs to reach 70 %. The network and label starts are within a few steps of the region where it moves fast.

### 4.3 What the starts' orbitals are worth: Jastrow-only optimization

Here the orbital rotations are frozen at the start and only the 2 × (2·norb − 1) Jastrow parameters are optimized
(L-BFGS, converged). This is the setting of hardware experiments that keep the orbitals fixed.

| molecule | truncated | cdf_lin | label | pre4 | rl4n29f |
|:--|--:|--:|--:|--:|--:|
| C2H3N_rxn2858_P | 42.9 | 51.7 | 77.1 | 75.3 | 75.3 |
| C2H3N_rxn2858_TS | 56.9 | 64.1 | 78.1 | 77.0 | 77.0 |
| C2H4O_rxn0724_TS | 33.6 | 56.1 | 78.8 | 77.9 | 77.9 |
| C2H4O_rxn2507_TS | 38.4 | 62.4 | 73.0 | 76.8 | 76.8 |
| C3H4_rxn2391_TS | 29.9 | 50.1 | 77.0 | 75.4 | 75.4 |

Our two orbital sets (label and network) are equally good. Lin et al.'s are 15–45 points worse. The network's
one-call lead over the label therefore comes from its Jastrow factors, not its orbitals. Once the Jastrows are
optimized from the label's orbitals, that lead is gone.

### 4.4 Final energies

Best over the whole run; the budget is the same for every start of a molecule:

| molecule | zero | random0 | truncated | cdf_lin | label | pre4 | rl4L | rl4n29f |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| C2H2O2_rxn3279_TS |  |  |  |  | 78.7 |  |  | 83.4 |
| C2H3N_rxn2858_P | -1.7 | 83.1 | 83.2 | 83.7 | 85.1 | 85.1 | 85.1 | 85.1 |
| C2H3N_rxn2858_TS | -3.3 | 85.8 | 87.5 | 83.4 | 87.0 | 88.7 | 88.7 | 88.7 |
| C2H4O_rxn0724_P |  |  | 84.0 |  | 87.0 |  |  | 85.6 |
| C2H4O_rxn0724_TS | -4.3 | 82.3 | 82.6 | 84.0 | 84.6 | 86.7 | 86.6 | 86.6 |
| C2H4O_rxn2507_P |  |  | 85.3 |  | 85.2 |  |  | 86.4 |
| C2H4O_rxn2507_TS | -5.4 | 82.2 | 84.8 | 81.7 | 83.5 | 85.5 | 85.5 | 85.5 |
| C2HNO_rxn4439_TS |  | 83.5 | 84.4 | 83.6 | 85.1 | 87.1 |  | 87.1 |
| C2HNO_rxn4440_TS |  | 84.1 | 85.5 | 84.4 | 87.3 | 87.5 |  | 87.5 |
| C2HNO_rxn4441_TS |  | 81.9 | 80.3 | 78.9 | 83.0 | 84.0 |  | 84.0 |
| C2HNO_rxn4442_TS |  |  |  |  | 83.0 |  |  | 83.4 |
| C3H4_rxn2391_P |  |  | 82.9 |  | 84.9 |  |  | 85.4 |
| C3H4_rxn2391_R |  |  | 87.6 |  | 87.4 |  |  | 88.3 |
| C3H4_rxn2391_TS | -8.4 | 80.0 | 83.7 | 81.5 | 83.9 | 83.5 | 83.4 | 83.7 |

## 5. Recommendations for Avijit

1. **Replace the objective.** BCH2 has no minimum, and no start can change that.
   * For validation on small systems: the exact energy, a drop-in with his own driver (`exact_model.py`).
   * For our molecules: the GPU engine here, exact to norb 18. Beyond that, our MPS energy engines; gradients there
     would need new work.
   * A polynomial surrogate has to be variational *in the UCJ state*. The truncated-BCH family is not.
2. **Use complex orbital rotations and two layers**, i.e. LUCJ.
   * With real rotations the first-order correlation is zero, and the ansatz holds 52–71 % of the H4 correlation.
   * With LUCJ, our network output is a native start: no projection, 0.2 s per molecule.
3. **Start from our network, or from our label if the network is unavailable.**
   * With exact gradients, both reach ~80 % of the CCSD correlation within ~10 evaluations at norb 15–18.
   * Lin et al.'s starts and random starts need ~80–100.
   * The zero start never leaves Hartree–Fock.
   * The network start ends higher than the label on most molecules, more often at norb 17–18.

## 6. Caveats

* **Precision.** Gradients are complex64/float32: ~1e-6 relative, 1e-6 Ha in the energy. No run reaches a gradient
  norm below ~1e-3, so "final" means the best within the budget, not a converged minimum.
* **Seeds.** One optimizer run per (molecule, start). Random starts: one seed, two for C2H3N. Budgets differ between
  molecules (§4).
* **The final orbital rotation stays the t1 rotation for all starts.** Optimizing it as well (Lin et al. do) adds
  norb² parameters per molecule.
* **Molecules.** 16 molecules: 9 at norb 15–16 (small_val), 4 at norb 17 (C2HNO TS) + 1 norb-17 P, 1 at norb 18
  (C2H2O2 TS). Small-sample statements are win counts, not significance tests.
* **Small systems.**
  * His seeds and starts are all-electron; the TS runs use our frozen-core active space (`active_mf`). He uses
    PySCF's CISD space, all electrons.
  * C2 under his driver stops at the 100-cycle cap.

## 7. Files

| file | what |
|:--|:--|
| `ucj_init/common.py` | loads his module (the folder's `__init__` imports a nonexistent `pyscf.ucj`); `ao2mo` shim and `active_mf` (his optimizer on our frozen-core active spaces); ffsim energy of his ansatz |
| `ucj_init/bch2_reproduce.py` | claims 1–3 + Taylor check + default runs, 7 small systems |
| `ucj_init/bch2_ts.py` | his optimizer on TS molecules from 5 starts; `projected_seed` (our LUCJ → his u0, g0) |
| `ucj_init/bch2_complex_check.py`, `rq2_test.py` | the two candidate repairs that fail |
| `ucj_init/exact_model.py`, `avijit_exact.py` | exact-energy `_LocalModel` for his kernel; runs |
| `ucj_init/ansatz_capacity.py`, `projection_loss.py` | ansatz capacity; projection loss |
| `ucj_init/lucj_grad_gpu.py`, `test_lucj_grad.py`, `check_ts_grad.py`, `profile_grad.py` | GPU adjoint LUCJ gradient, validation, profiling |
| `ucj_init/opt_exact.py` | runs: L-BFGS, trust region (`--method trn`), Jastrow-only (`--jastrow-only`) |
| `ucj_init/worker.sh`, `chain_worker.sh` | lock-based job queue |
| `ucj_init/report_exact.py`, `report_small.py`, `make_figs.py` | tables and figures |
| `rl_runs/ucj_init/exact_opt/*.jsonl` | every evaluation of every run |
| `pretrain/opt_true/results/ucj_init/exact_opt/*.json`, `*.npz` | summaries; start and best parameters |
