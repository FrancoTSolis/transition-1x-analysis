"""Avijit's UCJ-BCH2 optimizer (unmodified) on TS molecules, from several seeds.

The optimizer runs in our frozen-core active space (common.active_mf).  Seeds:
  default   : his CISD-singles kappa seed + relaxed-Hgg g seed (u0 = g0 = None)
  <cand>    : a project LUCJ candidate (label / rl4n29f / truncated / ...) projected onto his ansatz
              (real, one layer, spin-free):  R = polar(Re U_1)^T  (ffsim's U is R^T),  g = Z_1 on the square
              pattern (J_aa = J_ab = g in his spin-free Jastrow)
  random    : R = expm(A - A^T), A ~ 0.05 N(0,1);  g ~ 0.05 N(0,1) symmetric
Records E_BCH2 and |g| per cycle.  Every stationary point of E_BCH2 has E = E_HF (E_BCH2 - E_HF is exactly
quadratic in g at fixed U), so a seed can only decide how fast the run escapes.

    python -m ucj_init.bch2_ts --name C2H3N_rxn2858_TS --seed default --max-cycle 25 --out <json>
"""
from __future__ import annotations

import argparse
import json
import pickle
import time
from pathlib import Path

import numpy as np
import scipy.linalg

from ucj_init.common import ROOT, active_mf, load_avijit, make_ucj

CANDS = ROOT / "runs_ot" / "energy_tasks" / "task1_sqd_vs_lin" / "cands_smallval.pkl"


def projected_seed(U, Z):
    U1 = np.asarray(U)[0]
    W_, _, Vh = np.linalg.svd(U1.real)
    R = (W_ @ Vh).T
    n = R.shape[0]
    mask = np.eye(n, dtype=bool)
    p = np.arange(n - 1)
    mask[p, p + 1] = mask[p + 1, p] = True
    g = np.asarray(Z)[0] * mask
    return R, 0.5 * (g + g.T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--max-cycle", type=int, default=25)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    mod = load_avijit()
    d = np.load(ROOT / "rhf_hamiltonians" / f"{args.name}.npz")
    nocc = int(d["nelec_a"])
    mf = active_mf(d["one_body"], d["two_body"], nocc, float(d["constant"]))
    e_hf, e_ccsd = float(d["e_hf"]), float(d["e_ccsd"])
    ucj = make_ucj(mod, mf)
    ucj.verbose = 0
    ucj.max_cycle = args.max_cycle
    n = ucj.nmo
    u0 = g0 = None
    if args.seed == "random":
        rng = np.random.default_rng(0)
        A = 0.05 * rng.normal(size=(n, n))
        u0 = scipy.linalg.expm(A - A.T)
        G = 0.05 * rng.normal(size=(n, n))
        g0 = 0.5 * (G + G.T)
    elif args.seed != "default":
        tasks = pickle.load(open(CANDS, "rb"))["tasks"]
        U, Z = [(U, Z) for k, nm, U, Z, _ in tasks if nm == args.name and k[1] == args.seed][0]
        u0, g0 = projected_seed(U, Z)
    rec = dict(name=args.name, seed=args.seed, e_hf=e_hf, e_ccsd=e_ccsd, norb=n, nocc=nocc)
    if u0 is not None:
        rec["e_bch2_start"] = float(ucj.energy(g0, u0))
        rec["g_norm_start"] = float(np.linalg.norm(g0))
    # per-cycle |g| needs the iterate: wrap _LocalModel to record it
    gnorms = []
    orig = mod._LocalModel.__init__

    def rec_init(self, parent, g, rotation):
        orig(self, parent, g, rotation)
        gnorms.append(float(np.linalg.norm(g)))
    mod._LocalModel.__init__ = rec_init
    t0 = time.time()
    e_corr, g, U = ucj.kernel(g0=g0, u0=u0)
    mod._LocalModel.__init__ = orig
    rec.update(seconds=time.time() - t0, converged=bool(ucj.converged), cycles=len(ucj.history),
               e_corr_final=float(e_corr), g_norm_final=float(np.linalg.norm(g)),
               ccsd_corr=e_ccsd - e_hf,
               history=[dict(cycle=h["cycle"], e_corr=h["energy"] - e_hf, trust=h["trust_radius"],
                             rho=h["ratio"]) for h in ucj.history],
               g_norm_trace=gnorms[::max(1, len(gnorms) // 200)])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(rec, open(args.out, "w"), indent=1)
    print(f"{args.name} {args.seed}: converged={ucj.converged} cycles={len(ucj.history)} "
          f"E_corr={e_corr:+.4f} (CCSD {e_ccsd - e_hf:+.4f}) |g|={np.linalg.norm(g):.2f} {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
