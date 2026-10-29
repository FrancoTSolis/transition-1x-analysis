"""How much correlation can each UCJ variant hold?  Multi-start exact-energy minimization (L-BFGS, autograd)
on Avijit's test systems (all electrons, STO-3G), no final orbital rotation:

  real1    Avijit's ansatz: real orthogonal U (n(n-1)/2), one spin-free Jastrow g (J_aa = J_ab = g, dense)
  cplx1    the same with a complex unitary U (n^2)
  lucj2sq  our ansatz: 2 layers, complex U_r, J_aa on (p, p+1), J_ab on (p, p)   (no t1 rotation here)
  ucj2     2 layers, complex U_r, dense J_aa and J_ab

Starts: 6 random (log U ~ 0.3 N, J ~ 0.3 N) + for the complex variants the ffsim truncated t2 init.
Reports the best energy as % of the FCI correlation energy.

    python -m ucj_init.ansatz_capacity --out pretrain/opt_true/results/ucj_init/ansatz_capacity.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.optimize
import torch

from ucj_init.common import ref_energies, small_mf
from ucj_init.exact_model import fci_for_mf

torch.set_default_dtype(torch.float64)


def antiherm(v, n, real):
    K = torch.zeros(n, n, dtype=torch.complex128)
    iu = torch.triu_indices(n, n, 1)
    nt = iu.shape[1]
    K[iu[0], iu[1]] = v[:nt].to(torch.complex128)
    K = K - K.T
    if not real:
        iu2 = torch.triu_indices(n, n, 0)
        im = torch.zeros(n, n, dtype=torch.complex128)
        im[iu2[0], iu2[1]] = 1j * v[nt:nt + iu2.shape[1]].to(torch.complex128)
        im = im + im.T - torch.diag(torch.diagonal(im))
        K = K + im
    return K


def sym(v, n, mask=None):
    S = torch.zeros(n, n)
    iu = torch.triu_indices(n, n, 0)
    S[iu[0], iu[1]] = v
    S = S + S.T - torch.diag(torch.diagonal(S))
    return S if mask is None else S * mask


class Ansatz:
    def __init__(self, kind, fci):
        self.kind, self.fci = kind, fci
        n = fci.n
        self.n = n
        self.reps = 1 if kind in ("real1", "cplx1") else 2
        self.real = kind == "real1"
        self.nk = n * (n - 1) // 2 if self.real else n * n
        self.nj = n * (n + 1) // 2
        sq = torch.zeros(n, n)
        for p in range(n - 1):
            sq[p, p + 1] = sq[p + 1, p] = 1
        self.maa = sq if kind == "lucj2sq" else None
        self.mab = torch.eye(n) if kind == "lucj2sq" else None
        self.spinfree = kind in ("real1", "cplx1")
        self.per = self.nk + (self.nj if self.spinfree else 2 * self.nj)
        self.npar = self.reps * self.per

    def unpack(self, x):
        out, i = [], 0
        for _ in range(self.reps):
            K = antiherm(x[i:i + self.nk], self.n, self.real)
            i += self.nk
            if self.spinfree:
                g = sym(x[i:i + self.nj], self.n)
                Jaa = Jab = g
                i += self.nj
            else:
                Jaa = sym(x[i:i + self.nj], self.n, self.maa)
                i += self.nj
                Jab = sym(x[i:i + self.nj], self.n, self.mab)
                i += self.nj
            out.append((K, Jaa, Jab))
        return out

    def state(self, x):
        f = self.fci
        psi = f.hf
        occ = f.occ
        for K, Jaa, Jab in self.unpack(x):
            A = torch.linalg.matrix_exp(f.gen(K))                    # Lam(U), U = expm(K)
            psi = A.conj().T @ psi @ A.conj()                        # U^dag
            same = 0.5 * ((occ @ Jaa) * occ).sum(1)
            ph = same[:, None] + same[None, :] + occ @ Jab @ occ.T
            psi = psi * torch.exp(1j * ph)
            psi = A @ psi @ A.T                                      # U
        return psi

    def energy(self, x):
        return self.fci.energy(self.state(x))


def truncated_x(mf, ans):
    """ffsim truncated t2 init in this parameterization (complex variants only)."""
    import ffsim
    from pyscf import cc
    mycc = cc.CCSD(mf)
    mycc.kernel()
    op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(mycc.t2, n_reps=ans.reps)
    n = ans.n
    xs = []
    iu1 = np.triu_indices(n, 1)
    iu0 = np.triu_indices(n, 0)
    for r in range(ans.reps):
        K = ffsim.linalg.util.antihermitian_from_parameters(
            ffsim.linalg.util.unitary_to_parameters(op.orbital_rotations[r]), n) if False else None
        import scipy.linalg
        K = scipy.linalg.logm(op.orbital_rotations[r])
        K = 0.5 * (K - K.conj().T)
        kp = np.concatenate([K[iu1].real, K[iu0].imag])
        Jaa, Jab = op.diag_coulomb_mats[r]
        if ans.spinfree:
            xs += [kp, Jab[iu0]]
        else:
            if ans.maa is not None:
                Jaa = Jaa * ans.maa.numpy()
                Jab = Jab * ans.mab.numpy()
            xs += [kp, Jaa[iu0], Jab[iu0]]
    return np.concatenate(xs)


def minimize(ans, x0, maxiter=3000):
    maxiter = int(__import__("os").environ.get("CAP_MAXITER", maxiter))
    def fg(x):
        xt = torch.tensor(x, requires_grad=True)
        E = ans.energy(xt)
        (g,) = torch.autograd.grad(E, xt)
        return float(E), g.numpy().astype(np.float64)
    res = scipy.optimize.minimize(fg, x0, jac=True, method="L-BFGS-B",
                                  options=dict(maxiter=maxiter, maxfun=2 * maxiter, gtol=1e-8, ftol=1e-15))
    return float(res.fun), int(res.nfev), float(np.abs(res.jac).max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=["H4_lin_1.5", "H4_lin_2.0", "H4_sq_1.5", "LiH_1.6",
                                                     "LiH_3.0", "C2_1.25"])
    ap.add_argument("--kinds", nargs="*", default=["real1", "cplx1", "lucj2sq", "ucj2"])
    ap.add_argument("--n-random", type=int, default=6)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    res = []
    for name in args.systems:
        mf = small_mf(name)
        refs = ref_energies(mf)
        fci = fci_for_mf(mf)
        e_hf = mf.e_tot
        corr = refs["e_fci"] - e_hf
        for kind in args.kinds:
            ans = Ansatz(kind, fci)
            starts = [("random%d" % s, 0.3 * np.random.default_rng(s).normal(size=ans.npar))
                      for s in range(args.n_random)]
            if kind != "real1":
                try:
                    starts.append(("truncated", truncated_x(mf, ans)))
                except Exception as e:  # noqa: BLE001
                    print("truncated init failed", kind, e)
            runs = []
            t0 = time.time()
            for tag, x0 in starts:
                E0 = float(ans.energy(torch.tensor(x0)))
                E, nfev, gi = minimize(ans, x0)
                runs.append(dict(start=tag, E0=E0, E=E, pct=100 * (e_hf - E) / (e_hf - refs["e_fci"]),
                                 pct0=100 * (e_hf - E0) / (e_hf - refs["e_fci"]), nfev=nfev, grad_inf=gi))
            best = min(runs, key=lambda r: r["E"])
            rec = dict(name=name, kind=kind, npar=ans.npar, fci_corr=corr, ccsd_corr=refs["e_ccsd"] - e_hf,
                       best_pct=best["pct"], best_start=best["start"], runs=runs, seconds=time.time() - t0)
            res.append(rec)
            print(f"{name:11s} {kind:8s} npar {ans.npar:4d}  best {best['pct']:6.1f}% of FCI corr ({best['start']})  "
                  f"all: {' '.join('%.1f' % r['pct'] for r in runs)}  ({rec['seconds']:.0f}s)", flush=True)
            Path(args.out).parent.mkdir(parents=True, exist_ok=True)
            json.dump(res, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
