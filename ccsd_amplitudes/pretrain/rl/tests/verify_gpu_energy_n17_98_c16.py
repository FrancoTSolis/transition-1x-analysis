#!/usr/bin/env python3
"""[verifier] norb 17 (9,8): complex128 engine vs complex64 engine on C3HN_rxn2586_P (label, all1_t1, all4_t4)
on a 12 GB card (c16 squeezed with an explicit max_mem_gb), plus the memory plan the engine makes for norb 18
(9,9)/(10,8) under a 40 GB budget (tables only, alloc_state=False)."""
from __future__ import annotations
import json, pickle, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402
from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402

out = dict(plans={}, n17_98={})
for name in ["C4H2_rxn9388_P", "C2H6O_rxn2392_P"]:
    eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex64, max_mem_gb=40.0, alloc_state=False)
    out["plans"][name] = dict(eng.mem_plan, nq=eng.nq, K=eng.K, dim=eng.dim)
    print(name, json.dumps(out["plans"][name]), flush=True)
    eng.release(); del eng; torch.cuda.empty_cache()
tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks/n17_pre.pkl", "rb"))["tasks"]
sel = [t for t in tasks if t[1] == "C3HN_rxn2586_P"]
res = {}
for dt, mm in ((torch.complex64, None), (torch.complex128, 12.1)):
    torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    eng = LUCJEnergyGPU.from_npz("C3HN_rxn2586_P", dtype=dt, max_mem_gb=mm)
    for (nm, cand), _, U, Z, t1 in sel:
        t0 = time.time()
        E = eng.energy(U, Z, t1=t1)
        res.setdefault(cand, {})[str(dt)[6:]] = E
        print(f"{nm} {cand} {str(dt)[6:]}: E {E:.10f} corr% {100 * eng.corr_frac(E):.5f}  {time.time() - t0:.1f}s  tile {eng.T}"
              f"  peak alloc {torch.cuda.max_memory_allocated() / 1e9:.2f} GB", flush=True)
    eng.release(); del eng
for cand, d in res.items():
    d["c8_minus_c16"] = d["complex64"] - d["complex128"]
    print(f"{cand}: c8 - c16 = {d['c8_minus_c16']:+.2e}", flush=True)
out["n17_98"] = res
json.dump(out, open(ROOT / "pretrain/rl/tests/results/verify_gpu_energy_n17_98_c16.json", "w"), indent=1)
print("done")
