# optimize=True LUCJ pretraining that fits *and* generalizes: chemistry frame + slot-frame learned optimizer

*Follow-up to `optimize_true_learnability_study.md` (Sep 2026), 2 Oct 2026. Code: `pretrain/opt_true/`
(`frames.py`, `varpro.py`, `slotnet.py`, `train_slot*.py`, `eval_slot*.py`, `frame_eval.py`, `entanglement.py`),
`pretrain/rl/grpo_slot.py`, `expanse/`. Logs: `runs_ot/`, results JSON: `pretrain/opt_true/results/`.*

## 0. Summary

| Question | Answer |
|:--|:--|
| Can the network even fit the optimize=True labels? | **Yes.** With dropout really off, enough epochs, visible input and wide enough output bounds, the old MO-frame network memorizes the literal labels: 16 / 64 / 256 molecules reach the labels' own residual (U within 0.003 / 0.014 / 0.040 of the label). It does not transfer (val 0.87–0.93): the labels are not a function of $t_2$ (bit-near-identical inputs carry labels $dU\approx0.6$ apart). |
| What made the old pipeline look unable to fit? | `pretrain/train.py` never set `attention_dropout` (stayed 0.1 with `--dropout 0`); `PretrainingModel._init_weights` silently overwrote the zero-init of the residual heads; 40 epochs (memorization needs hundreds); $t_2$ enters at $2\times10^{-3}$ scale; tanh bound 1 below label entries up to ~2. The first two are fixed in this commit. |
| What generalizes? | Changing the **frame**, not the network. The optimize=True orbitals are atom-centred hybrids in a real sector. In a deterministic geometry-defined hybrid frame, a frame-consistent learned optimizer (slot model) trained only on ffsim's objective reaches **label quality on unseen molecules**: one call 0.37 (old net 0.79), 4 calls 0.336 = label 0.336, 8 calls **0.320 < label** on the leak-free split. |
| Energies (30 held-out molecules, exact statevector) | see §5 |
| RL | GRPO with exact-energy reward on the slot model runs on Expanse (§6). |
| Tensor-network reward | MPS+SQD of the 58-qubit circuit in the MO basis: 40 min per sample at $\chi=32$ and only 32 % correlation energy; not usable as a reward. Entanglement in the chemistry frame: §7. |

## 1. Why the old targets could be memorized but not learned

Memorization ladder (`runs_ot/fitlad_S*`): old MO-token network, supervised on the literal optimize=True
labels (phase-invariant $U$ loss + $Z$ loss), dropout 0 (incl. attention), constant LR, $t_2\times431$ input,
kscale 3 / zscale 2, zero-init heads.

| train set | epochs | sup. loss | train residual | $U$ distance to label | val residual |
|--:|--:|--:|--:|--:|--:|
| 16 | 1500 | 0.003 | 0.330 | 0.003 | 0.927 |
| 64 | 500 | 0.042 | 0.347 | 0.014 | 0.929 |
| 256 | 250 | 0.073 | 0.388 | 0.040 | 0.891 |
| 1269 | 100 | 0.262 | 0.601 | 0.267 | 0.865 (still falling) |

Label residual ≈ 0.32. Plateau escape time grows with $N$ (S256 left the 0.92 plateau after ~80 epochs).

## 2. The chemistry frame (real sector)

* Every optimize=True column has relative occupied/virtual phase $\mp\pi/4$: $U_k=\Phi_k O_k$,
  $\Phi_0=\mathrm{diag}(1_{occ},e^{-i\pi/4}1_{vir})$, $\Phi_1=\bar\Phi_0$, $O_k$ real orthogonal; then $\mathrm{Im}\,\mathrm{rec}=0$ exactly.
* Frame $B_k$ (`frames.py`, ~10 ms/molecule, all 30 205 molecules): valence AOs projected on the active MOs and
  Löwdin-orthonormalized (STO-3G valence = active space; MO-sign covariant), $p$ orbitals along the principal axes;
  rep 0 = bond-directed hybrids (lone pairs fill the exact orthogonal complement: planar 3-bond atoms get $p_\pi$),
  rep 1 = plain OAOs; chain order = DFS over heavy atoms with every tree bond chain-adjacent
  ($h_{A\to B}$ next to $h_{B\to A}$, $h_{A\to H}$ next to $s_H$).
* $U_k=\Phi_k B_k e^{A_k}$, $A_k$ real antisymmetric; $Z$ eliminated exactly (variable projection, `varpro.py`).

n29 val141, median ($Z$ by VarPro everywhere; `frame_eval.py`):

| start | residual |
|:--|--:|
| canonical exact-DF init (what optimize=True starts from) | 0.807 |
| chemistry frame, $A=0$ (no optimization) | 0.501 |
| + Adam on same-atom generators | 0.362 |
| + bonded generators | 0.319 |
| + unrestricted | **0.290** |
| optimize=True labels | 0.334 |

## 3. Slot-frame learned optimizer

Edge transformer whose tokens are the column pairs $(p,q)$ of the current iterate $U_t$; all inputs are expressed in
that same frame (projections of $t_2$ and of the residual on slot-pair dyads, the gradient of ffsim's objective,
Gram overlaps, $Z^*$, chain distance) plus chemistry features (element, hybrid type and direction, bonded / same-atom,
distance); it emits a real antisymmetric step $U_{t+1}=U_te^{\Delta K}$, $Z$ re-solved by VarPro, recycled $T$ times.
Trained on ffsim's objective only (no labels); Danskin gradient.

n29 (1269 train / 141 val; **val89** = the 89 val molecules without a train twin — the other 52 reactants have
bit- or near-identical copies in train):

| model | calls | val141 | val89 | train141 |
|:--|--:|--:|--:|--:|
| old MO-frame amortized net | 1 | 0.787 | — | 0.775 |
| slot model, canonical frame | 1 | 0.562 | 0.625 | 0.544 |
| slot model, chem frame (n29, d192) | 1 | 0.358 | 0.388 | 0.333 |
| slot model, chem frame (all 30k molecules) | 1 | 0.367 | **0.377** | 0.360 |
| slot model, chem frame, $T=4$ | 4 | 0.331 | 0.336 | 0.306 |
| slot model, chem frame, $T=8$ | 8 | **0.320** | **0.320** | 0.299 |
| optimize=True labels (L-BFGS ≤500 it.) | — | 0.334 | 0.336 | 0.317 |

By category ($T=4$ model, val141): reactants 0.319 → 0.307 (label 0.320), TS 0.402 → 0.351 (0.350),
products 0.377 → 0.319 (0.323) for one call → four calls. Mayer bond orders instead of distance bonds do not
improve the TS frame (0.526 → 0.529).

## 4. Fixes in the original pretraining code

* `pretrain/train.py`: `attention_dropout` now follows `--dropout` (was hard-wired 0.1).
* `pretrain/model.py`: zero-init of the residual heads restored after `_init_weights`.
* `pretrain/rl/grpo.py`: reward workers also pin `RAYON_NUM_THREADS` (ffsim's own pool; an unpinned worker
  takes ~22 cores); `pretrain/rl/hamiltonian.py`: MO signs aligned to the dataset's MOs.

## 5. Exact energies (30 held-out molecules, norb 15–16)

% of the CCSD correlation energy recovered, $(E_{HF}-E)/(E_{HF}-E_{CCSD})$, exact statevector (ffsim), square
connectivity, $n_{reps}=2$, $\lambda=0.005$, $t_1$ final rotation. None of these molecules was used in pretraining.
`pretrain/opt_true/results/energy_small_table.md` (all rows), JSON per molecule next to it.

| candidate | network calls / optimizer steps | median | mean | min | $t_2$ residual |
|:--|:--|--:|--:|--:|--:|
| canonical exact-DF init (optimize=True's start) | — | 14.3 | 16.1 | 4.9 | 0.866 |
| chemistry frame, exact Z | 0 | 34.1 | 30.4 | −64.3 | 0.448 |
| slot model, one shot (all sizes) | 1 call | 55.7 | 52.9 | 18.6 | 0.324 |
| slot model, 4 recycles (all sizes) | 4 calls | 58.9 | 57.3 | 25.2 | 0.272 |
| one shot + 100 in-frame Adam steps | 1 + 100 | 61.8 | 59.5 | 36.5 | 0.243 |
| chemistry frame + 300 in-frame Adam steps | 300 | 62.2 | 61.4 | 35.6 | 0.241 |
| **optimize=True labels** (L-BFGS ≤500 it.) | ≤500 | **67.6** | **64.3** | 36.4 | 0.253 |

* Residual is a weak proxy at this level: at equal residual the chemistry-frame basin is 2–5 points below the
  canonical/optimize=True basin (frame+Adam beats the label on 11/30 molecules, median −2.4).
* Not the layer order: swapping the two LUCJ layers (same residual, different circuit) moves energies by <1 point
  for every candidate (`*_swap` rows).
* This gap is what the energy RL stage is for (§6).

## 6. Energy RL (GRPO)

(filled in from `rl_runs/grpo_slot_v1`)

## 7. Tensor-network reward

Costs measured this session (single-threaded workers; ffsim's own thread pool must be pinned with
`RAYON_NUM_THREADS`, an unpinned worker takes ~22 cores):

| evaluator | system | cost per energy | notes |
|:--|:--|--:|:--|
| exact statevector | norb 15 (30 qubits) | ~200 core-s (scai1), 3 GB | 24 s wall with all threads on an Expanse node |
| exact statevector | norb 16 (32 qubits) | ~600 core-s, 10–12 GB | memory-bound: ~20 concurrent per 256 GB node |
| MPS ($\chi$=32) + SQD, MO basis | n29 (58 qubits) | 405 s MPS + 171 s sampling + 1773 s SQD | 126 unique of 2000 shots, **32 %** corr. energy |
| MPS ($\chi$=16) + SQD, MO basis | n29 (58 qubits) | 37 s + 249 s + 3536 s | 271 unique, **51 %**: lower fidelity, *higher* SQD energy |

The SQD energy tracks sample diversity, not LUCJ-parameter quality, at affordable $\chi$ (same conclusion as the
Sep study at 30 qubits), so MPS+SQD in the MO basis is not a usable RL reward.

Does the chemistry frame make the state MPS-friendly? Schmidt spectra of the exact label LUCJ state of
C2H3N_rxn2857_P (30 qubits, interleaved $\alpha\beta$ ordering, `entanglement.py`):

| basis | max entropy | max $\chi$ (1e-3 norm) | max $\chi$ (1e-4 norm) | mid-cut $\chi$ (1e-3 / 1e-4) |
|:--|--:|--:|--:|--:|
| MO (energy order) | 0.64 | 119 | 456 | 106 / 423 |
| frame rep 0 (hybrids, chain order) | 3.98 | 164 | 186 | 69 / 186 |
| frame rep 1 (OAOs, chain order) | 3.98 | 235 | 404 | 69 / 186 |

In the local frame the Hartree–Fock reference itself is entangled (bonds cross the cuts), so the head of the spectrum
is flatter while the tail decays faster: no clear win at 3 heavy atoms; untested for larger molecules.
(The $\alpha|\beta$ cut is basis-invariant, a sanity check: $S=0.387$, $\chi$ 29 / 93 in every basis.)
