# Can optimize=True LUCJ parameters be learned losslessly? Canonical gauge, optimizer multimodality, and what to train on instead

*Follow-up to `eigenvector_learnability_study.md` / `slides/deck_v3.md` (17 Sep 2026). All numbers
are from this repo's data; code in `gauge_study/compressed_canonical.py`, `generate_compressed_targets.py`,
`gauge_study/exp7–exp12`, `pretrain/train.py --loss-mode compressed_{sup,recon}`, `pretrain/rl/`.
Logs: `gauge_study/exp*_*.log`, `train_comp_*.log`, `rl_runs/`.*

---

## 0. Questions and verdicts

| Question | Verdict | Where |
|:--|:--|:--|
| Q1. If the gauge of the *optimize=True* (compressed-DF) labels is fixed at the source, is $(U,J)$ directly learnable? | **No.** Phase alignment is not the obstacle: ffsim's compressed optimizer is *multimodal*; restarting from a gauge copy of the same init gives an unrelated $U$, a different $Z$ and a 25 % different invariant content at the same residual. Canonical labels: $\kappa$ ridge $R^2$ 0.23, 1-NN orbit ratio 0.72 (random = 1). | §2, §3 |
| Q1′. Is there a *lossless* learnable formulation of the optimize=True task? | **Not as a single-valued regression target** — the argmin is a set (every restart lands in a distinct minimum, exp10). Lossless has to mean "produce *a* point at the label's objective value", which is achievable by (i) an amortized optimizer trained on the objective (recovers the smooth part, 0.92 → 0.79), (ii) multiple-hypothesis (winner-takes-all) heads scored by the residual at inference, (iii) predict-then-refine. And the minima are *not* energy-equivalent (7–11 points of correlation energy apart, exp12), so the tie-breaker should be energy, not residual. | §4, §5 |
| Q2. Which regularization $\lambda$? | $\lambda=0$ is **below Hartree–Fock** in energy (‖Z‖ blows up 10–20×). Plateau at $\lambda\in[0.005,0.01]$; 0.01 marginally best and safest. | §6 |
| Q3. Energy RL with a tensor-network reward? | GRPO pipeline built and running; the reward is the bottleneck: exact statevector 8–30 core-min/energy at norb 15–16; a *faithful* MPS of the LUCJ circuit needs $\chi\gtrsim 90$ already at 30 qubits (χ=64 loses 83 % of the norm), TN+SQD at small χ measures SQD's configuration recovery, not the parameters. | §7 |

---

## 1. Canonical-gauge compressed labels

ffsim's `optimize=True` path (`_double_factorized_t2_compressed`) is: exact truncated DF init (nested `eigh`,
LAPACK-arbitrary per-column phases, arbitrary basis in the $n_{occ}-n_{virt}$-dimensional zero eigenspace of
the quadratures, arbitrary sign of $v_0$) → L-BFGS on $\kappa=\log U$ and the allowed $Z$ entries. The objective
is invariant under the gauge group, but the L-BFGS path in $\kappa$-coordinates is not (column phases act on
$\kappa$ as $\log(UD)$, not as an isometry), so the output inherits and scrambles the init's gauge.

`canonical_exact_init` fixes every freedom of the init before optimizing: $v_0$ sign (largest-modulus entry
positive), zero-mode basis (pivoted QR of the gauge-invariant null projector), per-column phases (largest-modulus
entry real positive), column order (ascending $w$). `compress_from_init` then runs ffsim's own objective and
regularizer (`df_tensors_to_params` / `df_tensors_from_params_jax`) from that init. Labels were generated for the
n29 group (1410 molecules, $n_{occ}=16$, $n_{virt}=13$) with square connectivity $\lambda\in\{0,0.005\}$ and
all-to-all $\lambda=0.005$ (500 iterations, ~80 s per molecule per config on one core).

Label facts (square, $\lambda=0.005$): residual 0.929 (init) → 0.357 (q25 0.326, q75 0.381); L-BFGS rarely
converges in 500 iterations (1 %); $\|\Delta\kappa\|_F/\sqrt{2n}$ median 1.26; 29 % of molecules have eigenphases
of $U_0^\dagger U_{\rm opt}$ beyond $0.9\pi$ (the matrix log of the residual rotation sits on its branch cut);
$\|Z_{\rm opt}\|$ 1.20 vs 0.19 for the masked init and 0.63 for the full exact DF.

## 2. The optimizer, not the phase (exp8)

For 8 molecules, from the canonical init: (ii) restart from a *phase-rotated copy* of the same init; (iii) restart
from the canonical init of $t_2$ + 1 % noise (which moves the canonical init itself by 0.046). Distances are
phase-aligned and relative to $\|U\|$; $dZ$ relative to $\|Z\|$; $d\hat t_2$ relative to $\|t_2\|$. Medians:

| variant | resid | ‖Z‖ | conv. | (ii) dU | dZ | d t̂₂ | (iii) dU | dZ | d t̂₂ |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| exact-DF init (masked) | 0.932 | 0.19 | – | 0 | 0 | 0 | 0.046 | – | – |
| λ=0, 500 it | 0.344 | 5.73 | no | 1.10 | 1.29 | 0.30 | 0.90 | 0.76 | 0.25 |
| λ=0, 4000 it | 0.325 | 13.96 | no | 1.17 | 1.25 | 0.28 | 1.08 | 0.95 | 0.26 |
| λ=0.005, 500 it | 0.366 | 1.16 | no | 1.09 | 0.68 | 0.25 | 1.03 | 0.62 | 0.24 |
| λ=0.005, converged (942 it) | 0.358 | 1.17 | yes | 1.11 | 0.69 | 0.25 | 1.07 | 0.64 | 0.24 |
| + proximal 1e-3 ‖κ−κ₀‖² | 0.515 | 1.13 | yes | 0.57 | 0.85 | 0.24 | 0.54 | 0.62 | 0.24 |
| + proximal 1e-2 | 0.786 | 0.68 | yes | 0.13 | 1.11 | 0.11 | 0.07 | 0.55 | 0.09 |

$dU\approx1$ is the distance between two unrelated unitaries. Convergence does not help. The only way to make
the label a function of $t_2$ (proximal term) removes the gain (0.79 vs 0.36).

## 3. Learnability of the canonical labels (exp7)

Square, $\lambda=0.005$, all 1410 n29 molecules (`exp7_square_reg0.005_full.log`; the $\lambda=0$ and all-to-all
configs give the same picture). Closest-decile/all-pairs distance ratio (input contrast 0.881; lower = smoother,
1 = random): $Z_{\rm init}$ 0.431 · $B=\lambda_0v_0v_0^T$ 0.798 · $\hat t_2^{\rm opt}$ 0.861 · $Z_{\rm opt}$ 0.915 ·
$U_{\rm init}$ 0.967/0.965 (raw/phase-aligned) · $U_{\rm opt}$ **0.980/0.981** · $\Delta\kappa$ 0.999.
Phase alignment changes nothing; the optimizer made $Z$ 2× less smooth than the exact-DF $Z$ it started from.

Ridge on t2-PCA256, test $R^2$ pooled / **molecule-specific** (per-entry train means removed — the pooled number
is inflated by between-entry means): $Z_{\rm init}$ 0.959 / **0.847** · $\lambda_0v_0$ 0.609 / 0.566 ·
$\hat t_2^{\rm opt}$ 0.548 / 0.460 · $Z_{\rm opt}$ 0.817 / **0.283** · $\Delta Z$ 0.788 / 0.226 ·
$\kappa_{\rm opt}$ / $\Delta\kappa$ 0.228 / **0.224**, 0.230 / 0.228. The phase-free encoding
$H=U\,{\rm diag}(0..n{-}1)U^\dagger$ scores 0.885 pooled but 0.275 molecule-specific, and `eigh` of the
predicted $H$ gives a $U$ no closer to $U_{\rm opt}$ (1.26) than the train-mean $H$ (1.28).
1-NN transfer of $U_{\rm opt}$: 0.684 (raw) and 0.678 (phase-orbit) of the random-pair distance — the same
verdict under both metrics, unlike the exact-DF case of deck_v3 where the metric was the culprit.
Stratified by outer gap: $\kappa$ $R^2$ 0.13 / 0.25 / 0.21 for gap <0.05 / 0.05–0.2 / >0.2 — no stratum is learnable.

Training the Edge Transformer on these labels (`--loss-mode compressed_sup`, residual heads, MSE on
$\Delta\kappa$ and $\Delta Z$; 1124 molecules, 40 epochs): the $\Delta\kappa$ loss does not move (0.0529 →
0.0528), the $\Delta Z$ loss halves, the phase-aligned distance to $U_{\rm opt}$ stays at the init's 0.975, and
the residual of the predicted $(U,Z)$ is *worse* than the init (1.045 vs 0.918): regressing a mean over
multimodal labels moves $U$ off the init into nothing useful.

## 4. Amortizing the objective instead of the argmin

`--loss-mode compressed_recon`: residual heads on the canonical init, $U=U_0e^{\Delta\kappa}$,
$Z=Z_0+\Delta Z$ (zero-init heads ⇒ the network starts *at* the exact-DF init), trained on ffsim's objective
(resid² + λ-regularizer, per-molecule normalized) with no labels. Gauge and the $Z{=}0$ trap of v2's
reconstruction loss disappear by construction.

| n29 val (141 molecules), median t2 residual | |
|:--|--:|
| canonical exact-DF init (n_reps=2, square) | 0.919 |
| network, K=1, 40 epochs | **0.787** |
| per-molecule L-BFGS + proximal 1e-2 (deterministic) | 0.786 |
| per-molecule L-BFGS, λ=0.005 (multimodal labels) | 0.357 |
| supervised regression of those labels (40 ep) | 1.045 |
| K=4 / K=8 winner-takes-all heads, 40 epochs | 0.781 / 0.774 |
| K=4 heads with random κ-kicks 0.08, 40 epochs | **0.769** |
| network, K=1, full 30k dataset, 40 epochs (init 0.923) | 0.807 |
| network output + 32 L-BFGS iterations (exp11) | 0.40 (from the init: 132 iterations) |
| network output + 500 L-BFGS iterations | 0.286 (from the init: 0.289) |

The network reproduces exactly the deterministic part of the compressed-DF improvement — the same number the
proximally-constrained optimizer converges to — and plateaus from epoch 20. The rest of the gain is made of
basin choices that are not a continuous function of $t_2$.

**Why a single-valued target cannot be lossless (exp10).** Ten restarts (random phases + κ-kick 0.1) per
molecule give ten distinct minima (pairwise $d\hat t_2$ ≈ 0.20 ‖t₂‖) with residuals within 0.03 of each other:
the argmin set is a near-flat valley, not a few basins. Any continuous selection can only follow the valley
until it has to jump.

**The minima are not energy-equivalent (exp12).** For two norb=15 molecules, 9 restarts each, exact statevector
energies: residuals 0.244–0.278, correlation energy recovered 65.6–72.5 % and 59.1–69.8 % (spreads of 7 and 11
points; init 19.5 % / 13.1 %); correlation between residual and energy inside the valley −0.17 and −0.70. So
(a) *any* valley point is a 3–5× improvement over the exact-DF init, and (b) the residual is a weak tie-breaker
among valley points; energy should pick.

**Continuation labels do not help (exp13).** Chaining the optimization through greedy nearest neighbours in
$t_2$ (warm-starting molecule $k$ from molecule $k-1$'s solution) gives slightly *better* residuals (0.315 vs 0.327
fresh) but no smoother label field ($dU$ 0.99 vs 1.22, $d\hat t_2$ 0.98 vs 0.97 between consecutive molecules):
the n29 group is too sparse in $t_2$-space (nearest neighbours differ by 1.14 $\|t_2\|$) for any local
continuity to be exploitable.

**Lossless-by-construction options, in order of cost:**

1. *Multiple-hypothesis amortized optimizer* (`--residual-hyps K [--residual-hyp-kick s]`): K residual heads,
   loss $0.95\min_k \mathcal L_k + 0.05\,\mathrm{mean}_k \mathcal L_k$ (winner-takes-all with a relaxation),
   inference picks the head with the lowest residual — the objective is free at inference, so the selection is
   exact. Each head learns one smooth branch; the union is a piecewise-continuous selection, which is what a
   continuous regressor cannot be. Result (40 epochs, n29 val): K=1 0.787, K=4 0.781, K=8 0.774, K=4 with random per-head κ-kicks
   (0.08) **0.769** — a modest, monotone gain in K and head diversity, still far from the labels' 0.33
   (`train_comp_recon_n29_K{4,8,4kick}.log`).
2. *Predict-then-refine* (`gauge_study/exp11_refine.py`, 8 val molecules): from the network output, L-BFGS
   reaches residual 0.6 / 0.5 / 0.45 / 0.4 in 7 / 14 / 22 / 32 iterations vs 23 / 38 / 53 / 132 from the exact-DF
   init, and ends at 0.286 vs 0.289 after 500 — the same quality as the labels for ~4× less per-molecule
   optimization. This is the lossless recipe available today.
3. *Energy-scored hypotheses / RL*: score the K candidates by energy (exp12 says this matters at the 10 % level),
   or fine-tune the policy on energy directly (§7).
4. *Learned optimizer* (iterative network conditioned on the current iterate) — the general answer to amortizing
   a multimodal argmin; not attempted here.

## 5. What "aligned phases" can and cannot buy

The v3 result (invariant sufficient statistic of the exact DF) is lossless because the exact-DF label *is* a
function of $t_2$ up to gauge. For optimize=True the label is a function of $(t_2,\ \text{init gauge},\ \text{optimizer
path})$; removing the gauge from the init (this work) removes one source of randomness and leaves the other two,
which dominate: exp8 shows the same-molecule variability under a gauge copy of the init is as large as under a
1 % change of $t_2$. Hence the canonical labels are useful as *evaluation references* (their residual is the
target quality) but not as regression targets.

## 6. Regularization chosen by energy (exp9)

Exact statevector LUCJ energies (n_reps=2, with the $t_1$ final rotation) for the smallest molecules, % of CCSD
correlation energy recovered (norb=15, 6 molecules; extension to norb=16 in `exp9_sweep.log`):

| square | init | λ=0 | 0.001 | 0.002 | 0.005 | 0.01 | 0.02 | 0.05 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| median % | 22.9 | −271 | 62.2 | 62.9 | 69.9 | **71.3** | 67.2 | 50.0 |
| min % | 13.1 | −284 | 37.4 | 46.5 | 56.7 | **66.2** | 48.8 | 39.6 |
| residual | 0.824 | 0.194 | 0.204 | 0.216 | 0.233 | 0.272 | 0.362 | 0.569 |
| ‖Z‖ | 0.22 | 4.42 | 1.02 | 0.97 | 0.88 | 0.79 | 0.66 | 0.49 |
| **all-to-all** median % | 42.9 | −517 | 86.0 | 87.9 | 89.7 | **91.1** | 72.2 | 60.1 |

norb=16 molecules (square; 5–7 molecules per λ, the remaining norb=16 jobs hung in the process pool and were killed): init 13.6 %, λ=0.005
66.2 %, λ=0.01 **73.6 %**, λ=0.02 71.0 % (paired on 5 molecules: 66.2 / 73.9 / 70.1) — same ordering.
Unregularized compression is far below Hartree–Fock (Trotter error from the ‖Z‖ blow-up): the residual is not
the metric. Plateau at λ = 0.005–0.01, 0.01 marginally better and with the best worst case; λ ≥ 0.02 over-regularizes.
Regularized compression triples the correlation energy of the exact n_reps=2 truncation.

## 7. Energy reward and GRPO

**Costs measured** (`pretrain/rl/bench_energy.py`, `gauge_study/sqd_timing.py`, `mps_fidelity_check*.py`):
exact statevector (ffsim) 5–8 core-minutes at norb=15 (dim 4·10⁷), ~4× at norb=16, infeasible beyond.
MPS (quimb `CircuitMPS`) of the LUCJ circuit at 30 qubits: χ=64 keeps 17 % of the norm and
$|\langle HF|\psi\rangle|^2=0.002$ (exact 0.968), χ=128 46 % / 0.068, uncapped (max bond reached 87, cutoff
1e-10) exact. Interleaved αβ ordering and `CircuitPermMPS` are worse (norm < 1 %). The intermediate
$U^\dagger|HF\rangle$ states are rotated Slater determinants whose Schmidt rank grows exponentially with the
electron count at the cut; at norb≈29 (58 qubits) a faithful MPS is out of reach on our CPUs. SQD on the χ=32–64
samples still returns 91–97 % "correlation" because configuration recovery repairs particle-number-violating
samples — the reward would measure SQD, not the parameters — and its value swings 38 % → 97 % with the SQD
settings (batches, iterations).

**GRPO** (`pretrain/rl/grpo.py`, engineering from Transition-State-Generation-Flow's `compute_rl_loss`):
Gaussian policy over the free residual parameters $(\Delta\kappa,\Delta Z)$ with fixed σ, closed-form log-ratio and
KL, group-normalized advantages, PPO clipping, reference policy refreshed every K steps, anchor option, JSON-lines
diagnostics (reward mean/std, ratio, clip-active fraction, adv sign fractions, KL, grad norm, timings); rewards
$(E_{HF}-E)/(E_{HF}-E_{CCSD})$ from a CPU pool with cached Hamiltonians. Policy = pretrained Edge-Transformer
backbone + residual heads (zero-init ⇒ exact-DF init; from the `compressed_recon` checkpoint ⇒ amortized
compressed DF). Smoke test passed end-to-end. First real run (6 norb=15 molecules, exact reward, G=8, B=2,
σ=0.02, 5 min per step): the amortized-DF policy already recovers 33.0 % (train) / 23.8 % (val) of the correlation
energy vs 23.9 % / 13.1 % for the exact-DF init before any RL step; over 40 steps (3.3 h) the deterministic val
energy rose monotonically 23.8 → 25.4 → 25.9 → 26.0 → 26.7 → 27.5 → 28.3 → 29.0 → 29.8 % (evaluated every 5
steps; exact-DF init 13.1 %), while the group-mean train reward moved 0.30 → 0.31 (noisy: 2 of 5 molecules per
step). Two bugs found and fixed on the way: BLAS threads in spawned reward workers (load 190 on 48
cores) and dropout in the update pass (a train-mode forward changed μ enough at σ=0.02 to send the Gaussian
ratio to 0). Results continue in `rl_runs/norb15_exact_v1/log.jsonl`.

**Compute.** One GRPO step with G=8, B=8 needs 64 energies ≈ 10–30 core-hours at norb 15–16; 200 steps ≈
2 000–6 000 core-hours. This is the case for the CPU allocation the user offered, together with reward engineering
(fixed-subspace / projected energies, smaller active spaces as a curriculum) before scaling to norb≈29.

## 8. Reproduction

```bash
cd ccsd_amplitudes
# canonical compressed labels (sharded; ~80 s/molecule/config/core)
.lucj_venv/bin/python3 generate_compressed_targets.py --names-file gauge_study/names_n29_16_13.txt \
    --configs square_reg0.005 --shard 0 --n-shards 60
.lucj_venv/bin/python3 -m gauge_study.exp8_optimizer_determinism --n-mols 8 --n-procs 40   # exp8
.lucj_venv/bin/python3 -m gauge_study.exp7_compressed_canonical --config square_reg0.005    # exp7
.lucj_venv/bin/python3 -m gauge_study.exp9_reg_energy_sweep --n-procs 28                    # exp9 (needs rhf_hamiltonians/)
.lucj_venv/bin/python3 -m gauge_study.exp10_basins ; ... exp12_valley_energy
pretrain/.train_venv/bin/python3 -m pretrain.train --loss-mode compressed_recon [--residual-hyps 4] \
    --names-file gauge_study/names_n29_16_13.txt --epochs 40 --batch-size 12 --lr 5e-4 ...
pretrain/.train_venv/bin/python3 -m pretrain.rl.hamiltonian --names-file gauge_study/names_small_norb_le16.txt
GPU=2 ./run_grpo_small.sh                                                                   # GRPO
```
