# Overnight plan, 29 Oct 2026: NN initialization for Avijit's UCJ optimizer

Task (Fang): make Avijit's simulator (avijit_simulator/UCJ_optimizer.py, BCH2 trust-region Newton)
converge better on our TS data. Compare our network init with random/naive inits and with Lin et al.'s
(Wan-Hsuan) truncated (optimize=False) and compressed (optimize=True) inits.

## A. BCH2 diagnosis (UCJ_context.md says: reproduce claims 1-3, stop-and-report if they hold)
- [x] H4/LiH: claims 1-3 hold exactly; BCH2 = 2nd-order Taylor of the exact energy (ffsim check); runaway.
- [ ] all 7 small systems -> pretrain/opt_true/results/ucj_init/bch2_reproduce.json (running, log rl_runs/ucj_init/)
- [ ] TS molecules (frozen-core active space via fake RHF, common.active_mf): BCH2 from several seeds
      (Avijit default, projected label / rl4n29f / truncated, random) -> divergence is seed-independent.
Theory: for real K the first-order term vanishes, E_BCH2 - E_HF = g^T M(U) g exactly; every stationary point has
E = E_HF; M indefinite -> unbounded. No initialization can fix this.

## B. Fix = exact (variational) energy, LUCJ ansatz, Avijit's trust-region Newton driver
- B1 small systems (H4, LiH, C2): dense torch autodiff, exact Hv; Avijit's real ansatz and LUCJ.
- B2 TS molecules norb 15-16 (small_val): GPU adjoint gradient (H psi on GPU + transition 1-RDMs + expm chain
  rule), inits zero / random / truncated / cdf_lin / label / pre4 / rl4L / rl4n29f
  (runs_ot/energy_tasks/task1_sqd_vs_lin/cands_smallval.pkl). Metrics: evaluations to reach thresholds,
  energy at fixed budgets, converged flag, final energy.

## Rules
GPUs scai2: never 0/1; 5-7 used by others at start. scai3-7: <=2 fts GPUs per machine.
