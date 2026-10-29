"""Reproduce the three claims of UCJ_context.md ("Verified behaviour of the BCH2 objective") on H4, LiH, C2,
and check the optimizer's ansatz against an independent exact evaluation (ffsim).

Per system:
  1. E(g=0, U) = E_RHF for random orthogonal U.
  2. E(t g, U) - E_RHF = t^2 (E(g, U) - E_RHF)   (random g, U; t = 0.5, 2, 10).
  3. the g-g Hessian M(U) (dense, from the optimizer's g_hessp) at U = I and at the CISD seed: lowest eigenvalue.
  4. the exact energy of the same ansatz (ffsim) agrees with BCH2 to O(g^2): [E_exact(t g) - E_RHF] / t^2 ->
     E_BCH2(g) - E_RHF as t -> 0.
  5. the default run (CISD seed, 50 cycles): E_corr and |g| per cycle, converged flag, exact energy of the end point.

    python -m ucj_init.bch2_reproduce --out pretrain/opt_true/results/ucj_init/bch2_reproduce.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.linalg

from ucj_init.common import SMALL, ffsim_energy, load_avijit, mo_integrals, ref_energies, small_mf


def rand_orth(n, rng, scale=0.3):
    A = rng.normal(size=(n, n)) * scale
    return scipy.linalg.expm(A - A.T)


def rand_sym(n, rng, scale=0.1):
    A = rng.normal(size=(n, n)) * scale
    return 0.5 * (A + A.T)


def dense_hgg(ucj, g, U):
    model = ucj_mod._LocalModel(ucj, g, U)
    I = np.eye(ucj.ng)
    M = np.stack([model.g_hessp(I[i]) for i in range(ucj.ng)], 1)
    M = 0.5 * (M + M.T)
    # remove the phase (N^2 gauge) direction, which g_hessp projects out
    B = ucj._phase_complement()
    return B.T @ M @ B


def run_system(name, max_cycle=50, seed=0):
    rng = np.random.default_rng(seed)
    mf = small_mf(name)
    refs = ref_energies(mf)
    h1, eri = mo_integrals(mf)
    const = mf.mol.energy_nuc()
    ucj = ucj_mod.UCJ(mf)
    ucj.verbose = 0
    n, nocc = ucj.nmo, ucj.nocc
    e_hf = mf.e_tot
    out = dict(name=name, basis=SMALL[name][1], nmo=n, nocc=nocc, nkappa=ucj.nkappa, ng=ucj.ng, **refs)

    # 1. g = 0 gives E_RHF for any U
    d1 = [ucj.energy(np.zeros((n, n)), rand_orth(n, rng)) - e_hf for _ in range(3)]
    out["claim1_max_abs_E(g=0,U)-E_RHF"] = float(np.max(np.abs(d1)))

    # 2. homogeneity in g at fixed U
    U = rand_orth(n, rng)
    g = rand_sym(n, rng)
    base = ucj.energy(g, U) - e_hf
    rel = []
    for t in (0.5, 2.0, 10.0):
        rel.append(abs((ucj.energy(t * g, U) - e_hf) - t * t * base) / max(abs(t * t * base), 1e-300))
    out["claim2_E(g,U)-E_RHF"] = float(base)
    out["claim2_max_rel_dev_from_t^2"] = float(np.max(rel))

    # 3. lowest curvature of the (exactly quadratic) g-landscape
    seedU = ucj._initial_rotation()
    for tag, UU in (("I", np.eye(n)), ("cisd_seed", seedU), ("random", U)):
        w = np.linalg.eigvalsh(dense_hgg(ucj, np.zeros((n, n)), UU))
        out[f"claim3_Hgg_min_eig_U={tag}"] = float(w[0])
        out[f"claim3_Hgg_n_neg_U={tag}"] = int((w < -1e-10).sum())
        out[f"claim3_Hgg_max_eig_U={tag}"] = float(w[-1])

    # 4. ffsim exact energy of the same ansatz vs BCH2 (mapping check: BCH2 = 2nd-order Taylor of exact)
    chk = []
    for t in (1e-1, 3e-2, 1e-2):
        ex = ffsim_energy(h1, eri, const, nocc, t * g, U) - e_hf
        chk.append(dict(t=t, exact_over_t2=ex / t ** 2, bch2_over_t2=base))
    exm = ffsim_energy(h1, eri, const, nocc, 1e-2 * g, U, sign=-1) - e_hf
    out["claim4_exact_vs_bch2"] = chk
    out["claim4_sign_flip_same"] = bool(abs(exm / 1e-4 - chk[-1]["exact_over_t2"]) < 1e-6 * max(1, abs(base)))

    # 5. default run
    ucj.max_cycle = max_cycle
    hist = []
    t0 = time.time()
    e_corr, gfin, Ufin = ucj.kernel()
    out["run_seconds"] = time.time() - t0
    for h in ucj.history:
        hist.append(dict(cycle=h["cycle"], e_corr=h["energy"] - e_hf, trust=h["trust_radius"]))
    out["run_converged"] = bool(ucj.converged)
    out["run_cycles"] = len(ucj.history)
    out["run_e_corr_final"] = float(e_corr)
    out["run_g_norm_final"] = float(np.linalg.norm(gfin))
    grad, _ = ucj.gen_g_hop(gfin, Ufin)
    out["run_grad_inf_final"] = float(np.linalg.norm(grad, np.inf))
    out["run_exact_energy_at_end_minus_hf"] = float(ffsim_energy(h1, eri, const, nocc, gfin, Ufin) - e_hf)
    out["fci_corr"] = refs["e_fci"] - e_hf
    out["ccsd_corr"] = refs["e_ccsd"] - e_hf
    # |g| along the trajectory: re-derive from the stored history is not possible, so track e_corr only
    out["run_history"] = hist
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=list(SMALL))
    ap.add_argument("--max-cycle", type=int, default=50)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    res = []
    for s in args.systems:
        t0 = time.time()
        r = run_system(s, args.max_cycle)
        r["seconds"] = time.time() - t0
        res.append(r)
        print(f"{s:14s} nmo={r['nmo']:2d} E_corr FCI={r['fci_corr']:+.5f} CCSD={r['ccsd_corr']:+.5f} | "
              f"c1 {r['claim1_max_abs_E(g=0,U)-E_RHF']:.1e} c2 {r['claim2_max_rel_dev_from_t^2']:.1e} | "
              f"Hgg min eig I/seed {r['claim3_Hgg_min_eig_U=I']:+.3e}/{r['claim3_Hgg_min_eig_U=cisd_seed']:+.3e} "
              f"(neg {r['claim3_Hgg_n_neg_U=cisd_seed']}) | exact/t2 {r['claim4_exact_vs_bch2'][-1]['exact_over_t2']:+.6f} "
              f"bch2 {r['claim2_E(g,U)-E_RHF']:+.6f} | run: conv={r['run_converged']} cyc={r['run_cycles']} "
              f"E_corr={r['run_e_corr_final']:+.3f} |g|={r['run_g_norm_final']:.2f} "
              f"exact(end)-HF={r['run_exact_energy_at_end_minus_hf']:+.3f} ({r['seconds']:.0f}s)", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=1)


ucj_mod = load_avijit()

if __name__ == "__main__":
    main()
