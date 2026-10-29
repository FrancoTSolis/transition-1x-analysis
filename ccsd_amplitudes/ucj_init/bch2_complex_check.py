"""Does a complex orbital rotation make the BCH2 objective bounded?  For the one-layer spin-free UCJ with a
COMPLEX unitary U (exact energy, autograd): b = dE/dg and M = 1/2 d2E/dg2 at g = 0 (gauge direction removed).
BCH2 = E0 + b.g + g^T M g exactly (the second-order Taylor polynomial); it is bounded below iff M >= 0.
Real U: b = 0 identically.

    python -m ucj_init.bch2_complex_check
"""
import numpy as np
import scipy.linalg
import torch

from ucj_init.ansatz_capacity import Ansatz
from ucj_init.common import small_mf, ref_energies
from ucj_init.exact_model import fci_for_mf

torch.set_default_dtype(torch.float64)
for name in ["H4_lin_1.5", "H4_lin_2.0", "LiH_1.6", "C2_1.25"]:
    mf = small_mf(name)
    refs = ref_energies(mf)
    fci = fci_for_mf(mf)
    e_hf = mf.e_tot
    for kind in ["real1", "cplx1"]:
        ans = Ansatz(kind, fci)
        n = ans.n
        nj = ans.nj
        ones = np.zeros(nj)
        iu = np.triu_indices(n)
        ones[:] = 1.0
        P = scipy.linalg.null_space(ones[None] / np.linalg.norm(ones))   # remove the N^2 gauge
        rows = []
        for s in range(4):
            rng = np.random.default_rng(100 + s)
            k = 0.4 * rng.normal(size=ans.nk)
            def E_of_g(gv):
                return ans.energy(torch.cat([torch.as_tensor(k), gv]))
            g0 = torch.zeros(nj, requires_grad=True)
            E0 = E_of_g(g0)
            (b,) = torch.autograd.grad(E0, g0, create_graph=True)
            H = torch.stack([torch.autograd.grad(b[i], g0, retain_graph=True)[0] for i in range(nj)]).detach().numpy()
            b = b.detach().numpy()
            M = 0.5 * P.T @ H @ P
            bb = P.T @ b
            w = np.linalg.eigvalsh(0.5 * (M + M.T))
            bch2_min = float(E0) - 0.25 * bb @ np.linalg.solve(M, bb) if w[0] > 0 else -np.inf
            rows.append((np.linalg.norm(b), w[0], w[-1], int((w < -1e-10).sum()), len(w), bch2_min - e_hf))
        print(f"{name:11s} {kind}: |b| " + " ".join(f"{r[0]:.1e}" for r in rows) + " | min eig M " +
              " ".join(f"{r[1]:+.1e}" for r in rows) + f" | n_neg " + " ".join(f"{r[3]}/{r[4]}" for r in rows) +
              f" | BCH2 min - HF " + " ".join(f"{r[5]:+.3f}" for r in rows) +
              f" | FCI-HF {refs['e_fci'] - e_hf:+.3f}", flush=True)
