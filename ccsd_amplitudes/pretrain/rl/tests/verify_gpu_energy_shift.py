#!/usr/bin/env python3
"""[verifier] Robustness of the engine's float64 energy evaluator: the HF shift is exact algebra, so the energy
of a fixed state must not depend on the shift reference.  Rebuild the per-row coefficient table C and E0 at
runtime (engine object only; the deliverable is not modified) for shift occupations s_p in {HF (default),
0 (no shift, like an unshifted float64 evaluator such as ffsim's), 0.5}, and evaluate the same complex128 state.
"""
from __future__ import annotations
import argparse, json, pickle, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402
from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402


def rebuild(eng, h1, eri, shift):
    n = eng.norb
    hp = h1 - 0.5 * np.einsum("prrq->pq", eri)
    iu = np.triu_indices(n)
    Wm = eri[iu[0], iu[1]][:, iu[0], iu[1]]
    Wm = 0.5 * (Wm + Wm.T)
    lam, vec = np.linalg.eigh(Wm)
    thr = 1e-12 * np.abs(lam).max()
    pos, neg = lam > thr, lam < -thr
    V = np.concatenate([(vec[:, pos] * np.sqrt(lam[pos])).T, (vec[:, neg] * np.sqrt(-lam[neg])).T], 0)
    s = np.concatenate([np.ones(int(pos.sum())), -np.ones(int(neg.sum()))])
    nP = np.zeros(len(iu[0]))
    diag = iu[0] == iu[1]
    nP[diag] = 2.0 * shift[iu[0][diag]]
    v = V @ nP
    g = hp[iu] + (s * v) @ V
    E0 = float(hp[iu] @ nP + 0.5 * np.sum(s * v * v))
    Vext = torch.as_tensor(np.concatenate([V, g[None]], 0), device=eng.device, dtype=torch.float64)
    G = eng.g
    Vd = Vext[:, G["diagP"]]
    sh = torch.as_tensor(shift, device=eng.device, dtype=torch.float64)
    C = torch.empty_like(eng.C)
    step = 512
    for r0 in range(0, eng.dim, step):
        r1 = min(eng.dim, r0 + step)
        blk = Vext[:, G["exc_p"][r0:r1]].permute(1, 0, 2) * G["exc_s"][r0:r1, None, :]
        C[r0:r1, :, :-1] = blk.to(C.dtype)
        C[r0:r1, :, -1] = ((G["occf"][r0:r1] - sh) @ Vd.T).to(C.dtype)
    eng.C, eng.E0, eng.nq, eng.npos, eng.NLp = C, E0, V.shape[0], int(pos.sum()), V.shape[0] + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", default="baselines"); ap.add_argument("--name", default="C2N2_rxn3923_R")
    ap.add_argument("--cand", default="label"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{a.pkl}.pkl", "rb"))["tasks"]
    (_, _, U, Z, t1), = [t for t in tasks if t[0] == (a.name, a.cand)]
    d = np.load(ROOT / "rhf_hamiltonians" / f"{a.name}.npz")
    h1, eri = d["one_body"], d["two_body"]
    out = {}
    for dt in (torch.complex128, torch.complex64):
        eng = LUCJEnergyGPU.from_npz(a.name, dtype=dt, max_mem_gb=6.0)
        E_default = eng.energy(U, Z, t1=t1)
        psi = eng.state(U, Z, t1=t1).clone()
        n, k = eng.norb, eng.k
        res = dict(default=E_default)
        for tag, shift in [("hf_rebuilt", np.r_[np.ones(k), np.zeros(n - k)]), ("noshift", np.zeros(n)),
                           ("half", 0.5 * np.ones(n))]:
            rebuild(eng, h1, eri, shift)
            eng.psi.copy_(psi)
            res[tag] = float(eng._energy_of_psi())
        tag = str(dt)[6:]
        out[tag] = res
        print(tag, "  ".join(f"{k_}: {v:.12f} ({v - E_default:+.2e})" for k_, v in res.items()), flush=True)
        eng.release(); del eng; torch.cuda.empty_cache()
    json.dump(out, open(a.out, "w"), indent=1)
    print("->", a.out)


if __name__ == "__main__":
    main()
