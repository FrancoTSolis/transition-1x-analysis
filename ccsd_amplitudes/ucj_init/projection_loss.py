"""How much of an LUCJ start survives the projection onto Avijit's ansatz (real U, one layer, spin-free g)?
Exact energies (pretrain.rl.gpu_energy, all-to-all connectivity = his spin-free Jastrow J_aa = J_ab = g, no t1)
of the projected seeds (ucj_init.bch2_ts.projected_seed) vs the original LUCJ candidates.

    CUDA_VISIBLE_DEVICES=2 python -m ucj_init.projection_loss --out pretrain/opt_true/results/ucj_init/projection_loss.json
"""
import argparse, json, pickle
import numpy as np
from pretrain.rl.gpu_energy import LUCJEnergyGPU
from ucj_init.bch2_ts import CANDS, projected_seed

ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); args = ap.parse_args()
tasks = pickle.load(open(CANDS, "rb"))["tasks"]
names = sorted({n for _, n, *_ in tasks})
out = []
for name in names:
    eng = LUCJEnergyGPU.from_npz(name, max_mem_gb=4.0)
    for k, n, U, Z, t1 in tasks:
        if n != name or k[1] not in ("label", "rl4n29f", "truncated", "pre4"):
            continue
        E_orig = eng.energy(np.asarray(U), np.asarray(Z), t1=t1)
        R, g = projected_seed(U, Z)
        E_proj = eng.energy(R.T[None].astype(complex), g[None], t1=None, connectivity="all-to-all")
        # the same projection but keeping the complex orbital rotation of layer 1 (one layer, spin-free)
        E_l1c = eng.energy(np.asarray(U)[:1], g[None], t1=None, connectivity="all-to-all")
        rec = dict(name=name, cand=k[1], pct_orig=100 * eng.corr_frac(E_orig), pct_proj_real1=100 * eng.corr_frac(E_proj),
                   pct_layer1_complex=100 * eng.corr_frac(E_l1c))
        out.append(rec)
        print(f"{name:20s} {k[1]:9s} LUCJ {rec['pct_orig']:6.1f}%   projected (real, 1 layer) {rec['pct_proj_real1']:7.1f}%   "
              f"layer 1 complex {rec['pct_layer1_complex']:7.1f}%", flush=True)
    eng.release(); del eng
json.dump(out, open(args.out, "w"), indent=1)
