# Overnight plan, 29 Oct 2026: NN initialization for Avijit's UCJ optimizer

Task (Fang): make Avijit's simulator (avijit_simulator/UCJ_optimizer.py, BCH2 trust-region Newton)
converge better on our TS data. Compare our network init with random/naive inits and with Lin et al.'s
(Wan-Hsuan) truncated (optimize=False) and compressed (optimize=True) inits.
Report: docs/ucj_init/report.md, generated from docs/ucj_init/report_template.md by `python -m ucj_init.fill_report`
(also prints the statistics the Answer section quotes; check the text against them after any new run).

## Status (29 Oct, 07:00): done
- [x] BCH2 diagnosis: 7 small systems + C2H3N TS/P from 5 seeds. Unbounded; complex U and RQ2 do not fix it.
- [x] mapping to ffsim verified; projection loss: projected seeds land above HF.
- [x] avijit_exact.py (his kernel + exact model), ansatz_capacity.py (+ ansatz_capacity_c2.json): all systems.
- [x] C2 FCI reference corrected (pyscf nroots=1 found an excited state): `python -m ucj_init.fix_c2_fci` (idempotent)
      was run on all small-system JSONs; common.ref_energies now solves 4 roots.
- [x] lucj_grad_gpu.py validated; exact_opt: 125 runs (L-BFGS 15 molecules norb 15-18, trn 4, Jastrow-only 5).
- [x] report tables + figures regenerated (make_figs, fill_report); Answer numbers checked against fill_report stats.
- [x] commit + push; HANDOFF.md section 7; memory note.

Possible next steps: trn at norb 17-18; second random seed everywhere; optimize the final orbital rotation too;
an MPS-based gradient beyond norb 18.

## Rules
GPUs scai2: never 0/1. scai3-7: <=2 fts GPUs per machine (currently 2 on scai3, 2 on scai7).
Workers: ucj_init/worker.sh <jobs> <gpu>; launch remotely with `setsid nohup ... & disown` (not `&&`).
