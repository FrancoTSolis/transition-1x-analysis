"""A bounded, polynomial-cost completion of BCH2 ("RQ2") for Avijit's ansatz, tested against the exact energy.

His BCH2 energy is  E0 + <x0|H|x0> - <h0|xx0>  with x0 = X|0>, xx0 = P X x0 (P: reference + CISD space).
This is the second-order part of the Rayleigh quotient of the truncated state
    psi2 = |0> + i x0 - 1/2 xx0     (all in the CISD space),
    RQ2  = [E0 + <x0|H|x0> - <h0|xx0> + 1/4 <xx0|H|xx0>] / [1 + 1/4 <xx0|xx0>],
which is variational in the CISD space (RQ2 >= E_CISD), equals BCH2 + O(g^4), and needs only one more H
application on a vector his code already builds (xx0).  Here everything is evaluated in the dense FCI basis with
P = mask of determinants with total excitation level <= 2 (autograd, L-BFGS), and the optimum is compared with
the exact energy of the same parameters and with the exact-energy optimum.

    python -m ucj_init.rq2_test --out pretrain/opt_true/results/ucj_init/rq2_test.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.optimize
import torch

from ucj_init.ansatz_capacity import Ansatz
from ucj_init.common import ref_energies, small_mf
from ucj_init.exact_model import fci_for_mf

torch.set_default_dtype(torch.float64)


def sd_mask(fci):
    occ = fci.occ_np
    hf = np.zeros(fci.n)
    hf[:fci.k] = 1
    exc = ((occ - hf) > 0).sum(1)                 # excitation level of each alpha string
    tot = exc[:, None] + exc[None, :]
    return torch.as_tensor((tot <= 2).astype(float)).to(torch.complex128)


def x_apply(fci, A, g, v):
    """X v,  X = Lam(U) G Lam(U)^dag (G diagonal in the rotated basis with eigenvalue phase(g))."""
    w = A.conj().T @ v @ A.conj()
    w = w * fci.phase(g).to(torch.complex128)
    return A @ w @ A.T


def rq2_energy(ans, x, P):
    f = ans.fci
    (K, Jaa, Jab), = ans.unpack(x)
    A = torch.linalg.matrix_exp(f.gen(K))
    g = Jab
    ref = f.hf
    x0 = P * x_apply(f, A, g, ref)
    xx0 = P * x_apply(f, A, g, x0)
    psi = ref + 1j * x0 - 0.5 * xx0
    Hpsi = P * f.sigma(psi)
    return (torch.vdot(psi.reshape(-1), Hpsi.reshape(-1)).real / torch.vdot(psi.reshape(-1), psi.reshape(-1)).real)


def bch2_energy(ans, x, P):
    f = ans.fci
    (K, Jaa, Jab), = ans.unpack(x)
    A = torch.linalg.matrix_exp(f.gen(K))
    ref = f.hf
    x0 = P * x_apply(f, A, Jab, ref)
    xx0 = P * x_apply(f, A, Jab, x0)
    E0 = torch.vdot(ref.reshape(-1), f.sigma(ref).reshape(-1)).real
    return (E0 + torch.vdot(x0.reshape(-1), f.sigma(x0).reshape(-1)).real
            - torch.vdot(f.sigma(ref).reshape(-1), xx0.reshape(-1)).real)


def minimize(fun, x0, maxiter=4000):
    def fg(x):
        xt = torch.tensor(x, requires_grad=True)
        E = fun(xt)
        (g,) = torch.autograd.grad(E, xt)
        return float(E), g.numpy()
    r = scipy.optimize.minimize(fg, x0, jac=True, method="L-BFGS-B",
                                options=dict(maxiter=maxiter, maxfun=2 * maxiter, gtol=1e-9, ftol=1e-15))
    return r.x, float(r.fun), int(r.nfev), float(np.abs(r.jac).max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=["H4_lin_1.5", "H4_lin_2.0", "H4_sq_1.5", "LiH_1.6",
                                                     "LiH_3.0", "C2_1.25"])
    ap.add_argument("--n-start", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = []
    for name in args.systems:
        mf = small_mf(name)
        refs = ref_energies(mf)
        e_hf = mf.e_tot
        fci = fci_for_mf(mf)
        ans = Ansatz("real1", fci)
        P = sd_mask(fci)
        pct = (lambda E: 100 * (e_hf - E) / (e_hf - refs["e_fci"]))
        # sanity: BCH2 here = Avijit's E_BCH2 (checked through homogeneity and against his code in bch2_reproduce)
        runs = []
        for s in range(args.n_start):
            x0 = 0.3 * np.random.default_rng(s).normal(size=ans.npar)
            xr, Er, nfev, gi = minimize(lambda x: rq2_energy(ans, x, P), x0)
            Eex_at_rq2 = float(ans.energy(torch.tensor(xr)))
            xe, Ee, _, _ = minimize(lambda x: ans.energy(x), xr)       # exact optimum reached from the RQ2 optimum
            xe0, Ee0, _, _ = minimize(lambda x: ans.energy(x), x0)     # exact optimum from the same start
            Eb = float(bch2_energy(ans, torch.tensor(xr), P))
            runs.append(dict(start=s, rq2=Er, exact_at_rq2=Eex_at_rq2, exact_from_rq2=Ee, exact_from_start=Ee0,
                             bch2_at_rq2=Eb, nfev=nfev, grad_inf=gi,
                             pct=dict(rq2=pct(Er), exact_at_rq2=pct(Eex_at_rq2), exact_from_rq2=pct(Ee),
                                      exact_from_start=pct(Ee0))))
            print(f"{name:11s} start {s}: RQ2 min {pct(Er):6.1f}%  exact at RQ2 point {pct(Eex_at_rq2):6.1f}%  "
                  f"exact opt from there {pct(Ee):6.1f}%  exact opt from same start {pct(Ee0):6.1f}%  "
                  f"(CISD {pct(refs['e_cisd']):.1f}%)", flush=True)
        out.append(dict(name=name, **refs, runs=runs, cisd_pct=pct(refs["e_cisd"])))
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        json.dump(out, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
