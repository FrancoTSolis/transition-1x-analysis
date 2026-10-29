# Overnight plan, 29 Oct 2026: NN initialization for Avijit's UCJ optimizer

Task (Fang): make Avijit's simulator (avijit_simulator/UCJ_optimizer.py, BCH2 trust-region Newton)
converge better on our TS data. Compare our network init with random/naive inits and with Lin et al.'s
(Wan-Hsuan) truncated (optimize=False) and compressed (optimize=True) inits.
Report draft: docs/ucj_init/report.md (FILL markers = tables still to paste).

## Status (01:55)
- [x] A. BCH2 diagnosis: 7 small systems (bch2_reproduce.json) + C2H3N TS/P from 5 seeds (bch2_ts/). Unbounded.
- [x] complex-U BCH2 still unbounded (bch2_complex_check.py, output in the report); RQ2 bounded surrogate fails.
- [x] mapping to ffsim verified; projection loss (projection_loss.json): projected seeds are above HF.
- [~] avijit_exact.py (his kernel + exact model): H4s, LiH done; LiH_3.0, C2 running (CPU).
- [~] ansatz_capacity.py: H4s, LiH done; LiH_3.0, C2 running.
- [x] lucj_grad_gpu.py validated (rel 1e-6 vs ffsim; norb 15 vs gpu_energy 1e-6 Ha).
- [~] exact_opt queue: rl_runs/ucj_init/exact_opt/jobs_n15.txt (scai2 GPUs 2,3,4), jobs_n16.txt (scai3 GPUs 2,3;
      scai7 GPUs 4,7). Includes norb-17 C2HNO TS, trn (Avijit driver) and Jastrow-only (lbfgsJ) jobs.
      Deferred (pre-locked) norb-16 R/P jobs: rl_runs/ucj_init/exact_opt/deferred_keys.txt (rmdir the lock to run).
- [x] chained: jobs_n18.txt (scai7 after n16, H200 only), jobs_n16rp.txt (scai2 + scai3 after their queues)
- [ ] fill report tables: python -m ucj_init.report_exact [--method trn|lbfgsJ]; python -m ucj_init.report_small;
      python -m ucj_init.make_figs
- [ ] commit + push; HANDOFF.md section; memory note.

## Rules
GPUs scai2: never 0/1. scai3-7: <=2 fts GPUs per machine (currently 2 on scai3, 2 on scai7).
Workers: ucj_init/worker.sh <jobs> <gpu>; launch remotely with `setsid nohup ... & disown` (not `&&`).
