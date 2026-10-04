"""Coordinator benchmark: GPU exact engine at norb 17-18 on a big GPU (one energy per shape)."""
import pickle, sys, time, json
import numpy as np, torch
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac
d = pickle.load(open("runs_ot/energy_tasks/largeval_base.pkl", "rb"))
idx = json.load(open("rhf_dataset/_index.json"))
seen = set()
for (key, name, U, Z, t1) in d["tasks"]:
    shp = tuple(idx[name])
    if shp in seen or key[1] != "label":
        continue
    seen.add(shp)
    t0 = time.time()
    eng = LUCJEnergyGPU.from_npz(name)
    t1_ = time.time()
    torch.cuda.reset_peak_memory_stats()
    W, _, Vh = np.linalg.svd(U); E = eng.energy(W @ Vh, Z, t1=t1)
    t2 = time.time()
    E2 = eng.energy(W @ Vh, Z, t1=t1)
    t3 = time.time()
    print(f"{name} shape {shp}: setup {t1_-t0:.1f}s first {t2-t1_:.1f}s second {t3-t2:.1f}s peak {torch.cuda.max_memory_allocated()/1e9:.1f} GB corr% {100*corr_frac(E, eng.e_hf, eng.e_ccsd):.3f}", flush=True)
    del eng; torch.cuda.empty_cache()
