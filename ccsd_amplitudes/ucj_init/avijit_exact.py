"""Avijit's optimizer (RUCJ.kernel, unmodified) with the exact energy model (exact_model.install) on his test
systems (H4, LiH, C2; STO-3G, all electrons), from several seeds, next to the original BCH2 model.

Seeds (his ansatz: real orbital rotation U, one spin-free Jastrow g):
  default    his CISD-singles kappa seed (|kappa| = 0.01) + relaxed-Hgg g seed (|g| = 0.05), from the BCH2 model
  zero       u0 = I, g0 = 0
  random<s>  u0 = expm(A - A^T), A ~ 0.1 N(0,1);  g0 ~ 0.1 N(0,1) symmetric
  truncated  ffsim LUCJ from CCSD t2 (optimize=False, n_reps=1, all-to-all), projected: R = polar(Re U)^T, g = J_ab
  compressed ffsim optimize=True (compressed DF, n_reps=1, all-to-all), projected the same way

    python -m ucj_init.avijit_exact --out pretrain/opt_true/results/ucj_init/avijit_exact.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.linalg

from ucj_init.common import SMALL, load_avijit, mo_integrals, ref_energies, small_mf
from ucj_init.exact_model import fci_for_mf, install, ucj_state


def ffsim_seed(mf, optimize):
    import ffsim
    from pyscf import cc
    mycc = cc.CCSD(mf)
    mycc.conv_tol = 1e-9
    mycc.kernel()
    op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(mycc.t2, n_reps=1, optimize=optimize,
                                                    options=dict(maxiter=200) if optimize else None)
    U = op.orbital_rotations[0]
    W_, _, Vh = np.linalg.svd(U.real)
    R = (W_ @ Vh).T
    if np.linalg.det(R) < 0:
        R[:, -1] *= -1
    g = op.diag_coulomb_mats[0, 1]
    return R, 0.5 * (g + g.T)


def exact_energy_of(fci, g, U):
    import torch
    A0 = fci.lam(np.asarray(U).T)
    psi = ucj_state(fci, A0, torch.zeros(U.shape), torch.as_tensor(g))
    return float(fci.energy(psi))


def run(name, seeds, max_cycle):
    mf = small_mf(name)
    refs = ref_energies(mf)
    e_hf = mf.e_tot
    fci = fci_for_mf(mf)
    out = []
    for seed in seeds:
        mod = load_avijit()                         # fresh module: original BCH2 _LocalModel
        ucj = mod.UCJ(mf)
        ucj.verbose = 0
        n = ucj.nmo
        if seed == "default":
            u0 = ucj._initial_rotation()
            g0 = ucj.get_init_guess(u0)
        elif seed == "zero":
            u0, g0 = np.eye(n), np.zeros((n, n))
        elif seed.startswith("random"):
            rng = np.random.default_rng(int(seed[6:]))
            A = 0.1 * rng.normal(size=(n, n))
            u0 = scipy.linalg.expm(A - A.T)
            G = 0.1 * rng.normal(size=(n, n))
            g0 = 0.5 * (G + G.T)
        else:
            u0, g0 = ffsim_seed(mf, optimize=(seed == "compressed"))
        rec = dict(name=name, seed=seed, **refs, nmo=n, nocc=ucj.nocc)
        rec["E_exact_start"] = exact_energy_of(fci, g0, u0)
        rec["E_bch2_start"] = float(ucj.energy(g0, u0))
        # 1) original BCH2 model
        ucj.max_cycle = max_cycle
        t0 = time.time()
        e_corr, g, U = ucj.kernel(g0=g0.copy(), u0=u0.copy())
        rec["bch2"] = dict(converged=bool(ucj.converged), cycles=len(ucj.history), e_corr=float(e_corr),
                           g_norm=float(np.linalg.norm(g)), E_exact_end=exact_energy_of(fci, g, U),
                           seconds=time.time() - t0)
        # 2) exact energy model, same driver
        install(mod, fci)
        ucj2 = mod.UCJ(mf)
        ucj2.verbose = 0
        ucj2.max_cycle = max_cycle
        t0 = time.time()
        e_corr, g, U = ucj2.kernel(g0=g0.copy(), u0=u0.copy())
        grad, _ = ucj2.gen_g_hop(g, U)
        rec["exact"] = dict(converged=bool(ucj2.converged), cycles=len(ucj2.history), e_corr=float(e_corr),
                            E=float(ucj2.e_tot), g_norm=float(np.linalg.norm(g)),
                            grad_inf=float(np.abs(grad).max()), seconds=time.time() - t0,
                            trace=[h["energy"] - e_hf for h in ucj2.history])
        out.append(rec)
        ex = rec["exact"]
        print(f"{name:12s} {seed:10s} start exact {rec['E_exact_start'] - e_hf:+.5f} | BCH2: conv={rec['bch2']['converged']} "
              f"E_corr={rec['bch2']['e_corr']:+9.3f} (exact at end {rec['bch2']['E_exact_end'] - e_hf:+.4f}) | "
              f"EXACT: conv={ex['converged']} cyc={ex['cycles']:3d} E_corr={ex['e_corr']:+.5f} |grad|={ex['grad_inf']:.1e} "
              f"| FCI {refs['e_fci'] - e_hf:+.5f} CCSD {refs['e_ccsd'] - e_hf:+.5f}", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=["H4_lin_1.5", "H4_lin_2.0", "H4_sq_1.5", "LiH_1.6", "LiH_3.0",
                                                     "C2_1.25"])
    ap.add_argument("--seeds", nargs="*", default=["default", "zero", "random0", "random1", "truncated",
                                                   "compressed"])
    ap.add_argument("--max-cycle", type=int, default=100)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    res = []
    for s in args.systems:
        res += run(s, args.seeds, args.max_cycle)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        json.dump(res, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
