# Handoff: ML-LUCJ (ccsd_amplitudes/), state as of 28 Oct 2026

**Read in this order:**
1. this file;
2. `docs/followups_2026-10.md` — the latest results and the ranked next steps;
3. `slides/deck_v5.md` — the whole story in 26 slides;
4. `docs/chem_frame_pretraining.md` (how pretraining works), `docs/rl_larger_molecules.md` (RL at 15–29 orbitals);
5. `docs/optimize_true_learnability_study.md` (why the old label regression failed).

## 1. Where things stand

**The pipeline.** A network predicts LUCJ parameters (2 layers, square connectivity, spin-balanced; U and Z/J; the
final orbital rotation comes from t1) from the CCSD amplitudes of Transition1x molecules (STO-3G, frozen core,
15–44 orbitals, 30,205 molecules).

1. **Pretraining.** A slot model works in a chemistry frame: geometry-defined atom-centred hybrids. It is trained on
   ffsim's compressed-DF (optimize=True) objective itself, not on the labels (`runs_ot/slotall_T4/best.pt`,
   4 recycles).
2. **Fine-tuning.** GRPO on the LUCJ variational energy.

**Headline (variational energy, % of the CCSD correlation energy).** One policy (`rl_runs/grpo_n29_tn/policy_last.pt`,
tag `rl4n29f`) beats the per-molecule optimize=True labels at every size tested:

| set | energies | labels | policy | honest gain (3 seeds) |
|:--|:--|--:|--:|:--|
| norb 15–16 (9 molecules) | exact | 62.0 | 70.3 | |
| norb 17–18 (22) | exact | 61.7 | 68.6 | **+5.7 ± 0.7** (late-training mean; +7.4 was selection-inflated) |
| norb 19 (12, untouched) | exact | 66.7 | 68.9 | **about +2 ± 1** (≈ seed spread) |
| n29 test (24, untouched, 4 formulas) | MPS χ256 | 60.4 | 66.5 | 24/24 molecules better |

Labels run to convergence (5,000 iterations) leave the baseline unchanged (61.03 → 61.06 %).

**The caveat that matters most: SQD/QSCI.** Under Lin et al.'s noiseless protocol (arXiv:2511.22476; local checkout
`../wan-hsuan-lucj/`):
* A better variational energy does **not** give a better QSCI energy. The QSCI error follows the number of distinct
  sampled configurations (ρ −0.96).
* Our best QSCI policy is `rl4L`. It is level with Lin's compressed DF at 10⁵ samples (8.6 vs 8.0 mHa) and behind it at
  10⁶ (2.08 vs 0.69).
* Winning on SQD needs a QSCI-aware reward (§5).

**Other results:**
* **Reaction barriers.** 2-layer LUCJ overestimates barriers by ~15 kcal/mol against CCSD(T) (CCSD: 4.6).
* **Engines.** The dm TN engine is accurate to 44 orbitals. MPS sampling works at 58 qubits (χ512 ≡ exact for QSCI).
  The diagonalization (1–2 M determinants at norb 29) is the bottleneck.
* **NERSC.** The 2027 ERCAP GPU note (`docs/nersc_2027_gpu_budget.md`, model `docs/nersc_2027_cost_model.py`) asks
  for 10,000 node-hours. It went to Brandon; the deadline was 5 Oct 2026.

## 2. Code map

| path | what |
|:--|:--|
| `pretrain/opt_true/` | frames, slot models, `train_slot_all.py`, `eval_slot*.py`, `dftorch.py` (compressed-DF objective) |
| `pretrain/rl/grpo_slot.py` | GRPO driver. `--reward exact\|tn\|queue` (`--queue-kind exact\|tn`); `--init-policy`, `--prefix-T 3` (frozen recycles), `--tn-impl` |
| `pretrain/rl/policy_dump.py` | candidates → energy-task pickle. Tags `init:`, `init4:` (pretrained), policy checkpoints; `--labels-dir` adds the label |
| `pretrain/rl/gpu_energy.py` | exact LUCJ energy on one GPU, complex64 to ~1e-7 Ha. Up to norb 18 on 40 GB, norb 19 on ≥ 60 GB. `state()` returns the CI matrix |
| `pretrain/rl/tn_energy.py` | MPS energy, density-matrix engine (`method="dm"` on CUDA, `"zipup"` on CPU); verified |
| `pretrain/rl/tn_energy_v1.py` | frozen zip-up engine. **Every n29 number in the docs used it** (`--tn-impl v1`) |
| `pretrain/rl/tn_sampler.py` | perfect MPS sampling on the GPU, plus integrals in the samples' basis (F·S, not MO); tests in `tests/test_tn_sampler.py` |
| `pretrain/rl/reward_queue.py`, `launch_workers.sh` | shared-filesystem energy queue. GPU guard: ≤ 2 user GPUs per machine on scai3–7. `worker --kind exact\|tn --chi --tn-impl --zip-margin --block2-threads`; env `THREADS`, `WORKER_ARGS` |
| `pretrain/rl/queue_eval.py` | evaluate a task pickle through a queue (`--kind tn`, `--only-names`, `--skip-keys`) |
| `pretrain/rl/report_tables.py` | prints every table of `docs/rl_larger_molecules.md` |
| `pretrain/rl/hamiltonian.py`, `generate_compressed_targets.py` | active-space Hamiltonians (CPU, seconds); optimize=True labels (`--maxiter`, default 500) |
| `pretrain/rl/ccsd_t_refs.py`, `reaction_energies.py` | CCSD(T) references (cached JSON); reaction energies and barriers |
| `pretrain/followups/` | scripts of the 5 Oct follow-ups: `task1_sqd_lin.py` (QSCI harness, Lin protocol), `task2_*` (NOMAD/SPSA), `task3_*`, `task5_*`, `task6_*`, `task7_*`, `sampler_*` |
| `expanse/` | `expanse.py` (TOTP + ControlMaster helper: `run`/`put`/`get`), `sync.sh` (+ `sync_extra.txt`), slurm scripts |
| `slides/deck_v5.md` | Marp deck (+ `.html`/`.pdf`); figures from `slides/make_v5_figs.py` |

**Evaluation sets.**
* `pretrain/rl/small_val.txt` (9, norb 15–16), `large_val.txt` (22, norb 17–18), `gauge_study/names_norb19.txt` (12).
* `pretrain/rl/n29_val20.txt`, and `n29_test24.txt` (untouched; C3H5NO2 and C5H5N are in no RL set).
* Note: the norb-17 C2HNO reactant appears 4 times, and so does the norb-19 C2H3NO reactant.

**Candidate tags.**
* `label`, `labconv` (converged label);
* `pre1`/`pre4` (pretrained, 1 or 4 recycles);
* `rl1s`/`rl4s` (small-set RL), `rl4L` (`grpo_large_rl4` step 25), `rl4Lf` (step 50);
* `rl4n29`/`rl4n29f` (`grpo_n29_tn` best/last);
* `cdf_lin` (Lin's unregularized compressed DF).

**Results.**
* JSON in `pretrain/opt_true/results/` (`energy_*.json`, `followups/<task>/`); RL logs in `rl_runs/*/log.jsonl`
  and `rl_runs/followups/`.
* Gitignored: checkpoints (`*.pt`), labels (`rhf_targets_compressed*`), Hamiltonians, queues (`rl_queue/`), MPS
  scratch, sample `.npz` files and task pickles.

## 3. Compute rules (from Fang; keep)

* **scai2** (local, 8× TITAN Xp 12 GB, 48 cores): free to use, except **never GPUs 0/1**. GPU 0 runs Fang's
  `jevplan` server.
* **scai1** (8× V100 16 GB, 80 cores): free, but skip busy GPUs. It is often saturated by other users; check the load
  first.
* **scai3–scai7**: at most **2 GPUs per machine for user fts, counting all of Fang's jobs**.
  * scai4 and scai6 are usually already at 2 (his gpu4pyscf jobs), so they are off-limits.
  * scai7 H200s are shared with other users' 70–130 GB jobs: cap memory (`--max-mem-gb`).
  * scai3 ssh hung on 5 Oct; it was reachable again on 28 Oct.
* **CPU-heavy work** only on scai1/scai2, or on Expanse.
* **Expanse** (account cla361, shared):
  * keep jobs economical (the 5 Oct follow-ups used 378 core-hours);
  * use `expanse/expanse.py`;
  * never print or copy the TOTP seed;
  * never touch the remote `fermi_arc/` directory, the conda env `negf`, or the original batch_screen project.
* **Processes.** Never kill processes you did not start. As of 28 Oct nothing from this project runs on scai2 or scai7.
* **Marp.** The system node is v6 and cannot run the installed marp-cli 4.5, which needs Node ≥ 22.12. PDF export also
  needs a Chrome (set `CHROME_PATH`). Download both into a scratch directory; they are not installed.

## 4. Pitfalls already hit

* **Lazy imports.** TN workers import the engine on their first task. Swapping `tn_energy.py` mid-run silently mixes
  engine versions. Pin `--tn-impl` for any result set.
* **Heartbeats.** They come from a child process (block2 holds the GIL; a heartbeat thread starved). TN collectors use
  `stale_s=1800`. Old orphan results in `rl_queue/*/done` are harmless.
* **Backgrounding.** `cmd | tail && nohup x &` backgrounds the whole `&&` list; this once submitted an evaluation twice.
  A wait loop that runs `pgrep -f "<pattern>"` on its own command line matches itself and never exits.
* **GPU memory.** An exact norb-18 worker next to TN workers on a 40 GB A100 ran out of memory. Size workers by the
  nvidia-smi footprint, which is 2–3× torch's allocated memory.
* **Absolute values.** MPS χ256 (v1, zip margin 1.5) is 0.4–1.5 points below converged at n29. Comparisons between
  candidates were stable to < 0.1 point.
* **Claims.** Checkpoints picked on the reported set inflate gains by ~1.5 points. Run ≥ 3 seeds, and report a fixed
  step, a late-training mean, or a separate test set.
* **Session limit.** Six long parallel agents hit the account's session limit after ~3 h on 5 Oct. Their detached jobs
  kept running. Keep agents lean and write state to disk so they can be resumed.

## 5. Next steps (ranked; details in `docs/followups_2026-10.md` §3)

1. **Decide the target metric** (variational energy or QSCI). That choice fixes the RL reward.
2. **Count distinct sampled strings at n29** for `label`, `cdf_lin` and `rl4n29f` (`tn_sampler`, χ512, ≤ 20 s per 10⁶
   samples). At norb 15–17 this count predicted the QSCI ranking.
3. **Prototype a QSCI-aware reward at norb 15–17**: the energy plus a coverage term, or QSCI at 10⁵ samples
   (~10 s GPU + 20–50 s on 4 CPU threads per reward).
4. **Build and benchmark the GPU SQD solver (SBD).** It gates norb-29 QSCI and a QSCI reward beyond norb 17.
5. **Bring `docs/rl_larger_molecules.md` in line with the seed-averaged numbers** above (not done yet; HANDOFF and
   deck_v5 already use them).
6. **NERSC-relevant measurements on an A100** (scai3 is reachable again):
   * the anchor at 29 and 44 orbitals;
   * χ64 vs χ128 group ranking at 33–44 orbitals.
7. **Separate batch size from reward precision at n29**: a v1 χ64 run with 6 molecules per step, ~2 h on one Expanse
   node.
8. **Best-of-k labels at norb 19**, which settles the fragile norb-19 margin.
9. **Heavy-hex variant of the policy → IBM QPU runs**, following Lin et al.'s hardware protocol.

The table of 5 Oct tasks (all done, except SBD) is in `docs/followups_2026-10.md` §2; each task has a report in
`docs/followups/`.

## 6. Tree state (28 Oct)

* Everything above is committed and pushed: deck_v5 in `9203fa1`, and this handoff in the commit right after it.
* Untracked: `../orbital_ordering/t1x_formula_split.json` at the repository root (created 28 Oct, not by this
  project's sessions; left alone). Also `../NN_LUCJ_full_outline/` and `../Transition-State-Generation-Flow/`, which are
  other projects.
