import pickle, sys, time, json
import numpy as np, torch
sys.path.insert(0, "/xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes")
from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac
idx = json.load(open("rhf_dataset/_index.json"))
d = pickle.load(open("runs_ot/energy_tasks/n17_pre.pkl", "rb"))
seen = set()
for (key, name, U, Z, t1) in d["tasks"]:
    shp = tuple(idx[name])
    if shp in seen or key[1] != "label": continue
    seen.add(shp)
    torch.cuda.reset_peak_memory_stats()
    try:
        eng = LUCJEnergyGPU.from_npz(name)
        W, _, Vh = np.linalg.svd(U); t0 = time.time(); E = eng.energy(W @ Vh, Z, t1=t1); t1_ = time.time(); E = eng.energy(W @ Vh, Z, t1=t1); t2 = time.time()
        print(f"{name} {shp}: {t1_-t0:.1f}s / {t2-t1_:.1f}s peak {torch.cuda.max_memory_allocated()/1e9:.1f} GB corr% {100*corr_frac(E, eng.e_hf, eng.e_ccsd):.3f}", flush=True)
        del eng
    except Exception as e:
        print(f"{name} {shp}: FAILED {type(e).__name__}: {str(e)[:150]}", flush=True)
    torch.cuda.empty_cache()
