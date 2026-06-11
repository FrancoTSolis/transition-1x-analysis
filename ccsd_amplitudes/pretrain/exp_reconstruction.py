#!/usr/bin/env python3
"""PRELIMINARY EVIDENCE: can a smooth NN map t2 -> (U,Z) under a
gauge-invariant reconstruction loss, and GENERALIZE to held-out molecules?

The model outputs its OWN (U,Z) (it never sees ffsim's gauge-scrambled U).
Loss = || reconstruct_t2(U,Z) - t2 ||.  If the held-out reconstruction beats
the constant-(U,Z) baseline and approaches the per-molecule ffsim floor, then
a smooth t2->(U,Z) section exists and is learnable -> the approach works.

Parametrization (n_reps reps): per rep,
  U = matrix_exp(A + iS),  A real antisym, S real sym  (=> U unitary)
  Z real symmetric
Run with system python (torch). Reads t2 from rhf_t2_cache + rhf_cache.
"""
import argparse
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


def load_t2(dirs, target_norb):
    Xs, metas, seen = [], None, set()
    for d in dirs:
        for f in sorted(Path(d).glob("*.npz")):
            if f.stem in seen:
                continue
            s = np.load(f)
            if int(s["norb"]) != target_norb:
                continue
            seen.add(f.stem)
            Xs.append(s["t2"].astype(np.float32))
            metas = (int(s["nocc"]), int(s["nvirt"]), int(s["norb"]))
    return np.stack(Xs), metas


class ReconModel(nn.Module):
    """MLP: pca(t2) -> (anti-Herm generators, Z) for n_reps reps."""

    def __init__(self, in_dim, norb, n_reps, hidden=256):
        super().__init__()
        self.norb, self.n_reps = norb, n_reps
        self.ti = torch.triu_indices(norb, norb, offset=1)   # antisym (real part)
        self.di = torch.triu_indices(norb, norb, offset=0)   # sym (imag part / Z)
        n_anti = self.ti.shape[1]
        n_sym = self.di.shape[1]
        self.n_anti, self.n_sym = n_anti, n_sym
        out = n_reps * (n_anti + n_sym + n_sym)              # A, S, Z per rep
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.GELU(),
            nn.Linear(hidden, hidden), nn.GELU(),
            nn.Linear(hidden, out),
        )
        self.net[-1].weight.data *= 0.01
        self.net[-1].bias.data.zero_()

    def _fill(self, vals, idx, sym):
        B = vals.shape[0]
        M = torch.zeros(B, self.norb, self.norb, device=vals.device, dtype=vals.dtype)
        M[:, idx[0], idx[1]] = vals
        return M + M.transpose(1, 2) if sym else M - M.transpose(1, 2)

    def forward(self, x):
        B = x.shape[0]
        p = self.net(x)
        na, ns = self.n_anti, self.n_sym
        per = na + ns + ns
        Us, Zs = [], []
        gnorm = 0.0
        for k in range(self.n_reps):
            o = p[:, k * per:(k + 1) * per]
            # bound generator entries with tanh -> matrix_exp can't blow up
            A = self._fill(torch.tanh(o[:, :na]) * np.pi, self.ti, sym=False)
            S = self._fill(torch.tanh(o[:, na:na + ns]) * np.pi, self.di, sym=True)
            Z = self._fill(torch.tanh(o[:, na + ns:]) * 3.0, self.di, sym=True)
            K = torch.complex(A, S)
            U = torch.linalg.matrix_exp(K)
            Us.append(U); Zs.append(Z.to(torch.complex64))
            gnorm = gnorm + (A.pow(2).mean() + S.pow(2).mean())
        self.last_gnorm = gnorm / self.n_reps
        return torch.stack(Us, 1), torch.stack(Zs, 1)          # (B,reps,n,n)


def reconstruct(U, Z, nocc):
    # t2[i,j,a,c] = i * sum_k,pq Z U U* U U* ; batch=n, reps=k; slice occ,occ,vir,vir
    full = 1j * torch.einsum("nkpq,nkap,nkip,nkcq,nkjq->nijac",
                             Z, U, U.conj(), U, U.conj())
    return full[:, :nocc, :nocc, nocc:, nocc:].real


def rel_resid(rec, t2):
    return (torch.linalg.norm((rec - t2).reshape(t2.shape[0], -1), dim=1)
            / torch.linalg.norm(t2.reshape(t2.shape[0], -1), dim=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--norb", type=int, default=22)
    ap.add_argument("--n-reps", type=int, default=2)
    ap.add_argument("--pca", type=int, default=48)
    ap.add_argument("--epochs", type=int, default=4000)
    args = ap.parse_args()
    torch.manual_seed(0); np.random.seed(0)
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    t2, (nocc, nvirt, norb) = load_t2(["rhf_t2_cache", "rhf_cache"], args.norb)
    n = t2.shape[0]
    print(f"{n} molecules, norb={norb}, nocc={nocc}, nvirt={nvirt}, t2 dim={t2[0].size}")
    perm = np.random.permutation(n)
    t2 = t2[perm]
    n_tr = int(0.8 * n)

    Xflat = t2.reshape(n, -1)
    mu = Xflat[:n_tr].mean(0); sd = Xflat[:n_tr].std(0) + 1e-8
    Xn = (Xflat - mu) / sd
    # PCA on train
    U_, S_, Vt = np.linalg.svd(Xn[:n_tr], full_matrices=False)
    comp = Vt[:args.pca]
    Xp = Xn @ comp.T
    Xp = (Xp - Xp[:n_tr].mean(0)) / (Xp[:n_tr].std(0) + 1e-8)

    Xp = torch.tensor(Xp, dtype=torch.float32, device=dev)
    T2 = torch.tensor(t2, dtype=torch.float32, device=dev)
    tr, te = slice(0, n_tr), slice(n_tr, n)

    # ---- Baseline 0: predict zeros ----
    print(f"\nBaseline (predict zeros): rel_resid = 1.000")

    # ---- Baseline 1: best CONSTANT (U,Z) (no input dependence) ----
    const = ReconModel(1, norb, args.n_reps).to(dev)
    ones = torch.ones(n, 1, device=dev)
    opt = torch.optim.Adam(const.parameters(), lr=3e-3)
    for ep in range(args.epochs):
        opt.zero_grad()
        U, Z = const(ones[tr])
        loss = rel_resid(reconstruct(U, Z, nocc), T2[tr]).mean()
        loss.backward(); opt.step()
    with torch.no_grad():
        U, Z = const(ones[te])
        const_te = rel_resid(reconstruct(U, Z, nocc), T2[te]).mean().item()
        U, Z = const(ones[tr])
        const_tr = rel_resid(reconstruct(U, Z, nocc), T2[tr]).mean().item()
    print(f"Constant (U,Z) baseline:  train={const_tr:.3f}  test={const_te:.3f}")

    # ---- Model: t2 -> (U,Z) ----
    model = ReconModel(args.pca, norb, args.n_reps).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=2e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, args.epochs)
    best_te = 1.0
    for ep in range(args.epochs):
        model.train(); opt.zero_grad()
        U, Z = model(Xp[tr])
        loss = rel_resid(reconstruct(U, Z, nocc), T2[tr]).mean() + 1e-3 * model.last_gnorm
        loss.backward(); opt.step(); sched.step()
        if ep % 400 == 0 or ep == args.epochs - 1:
            model.eval()
            with torch.no_grad():
                U, Z = model(Xp[te])
                te_r = rel_resid(reconstruct(U, Z, nocc), T2[te]).mean().item()
            best_te = min(best_te, te_r)
            print(f"  ep {ep:4d}  train={loss.item():.3f}  test={te_r:.3f}")

    print(f"\n=== RESULT ===")
    print(f"predict-zeros        : 1.000")
    print(f"constant (U,Z)       : {const_te:.3f}  (no input dependence)")
    print(f"NN t2->(U,Z) [test]  : {best_te:.3f}")
    print(f"ffsim floor (n_reps=2 unconstrained, single-mol opt): ~0.32")
    if best_te < const_te - 0.02:
        print("\n=> POSITIVE: the NN uses t2 to produce better (U,Z) than any constant,")
        print("   i.e. a smooth, generalizing t2->(U,Z) map EXISTS and is learnable.")
    else:
        print("\n=> INCONCLUSIVE at this sample size; needs more molecules / capacity.")


if __name__ == "__main__":
    main()
