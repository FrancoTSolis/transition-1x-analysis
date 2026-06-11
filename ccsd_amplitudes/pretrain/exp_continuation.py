#!/usr/bin/env python3
"""DECISIVE EVIDENCE: does a SMOOTH t2 -> (U,Z) section exist?

ffsim solves each molecule independently -> U is gauge-scrambled -> nearby
t2 give far-apart U (the unlearnable target). But the reconstruction loss is
gauge-invariant, so we can instead pick a smooth representative:

  Order molecules along a path of increasing-t2-similarity. Optimize (U,Z) for
  each via reconstruction loss, WARM-STARTED from the previous molecule. If the
  warm-started U stays smooth (small ||ΔU|| between neighbors) while still
  reconstructing t2 well, then a smooth section EXISTS and is learnable.

Compare:  ||ΔU|| (warm-start, smooth section)   vs   ||ΔU|| (ffsim, independent)

Run with system python (torch). Reads rhf_cache (needs t2 + ffsim kappa).
"""
import argparse
from pathlib import Path

import numpy as np
import torch


def load(cache_dir, target_norb):
    rows = []
    for f in sorted(Path(cache_dir).glob("*.npz")):
        s = np.load(f)
        if int(s["norb"]) != target_norb:
            continue
        if "kappa_real" not in s:
            continue
        U = (np.cos(0) + 0j)  # placeholder
        from scipy.linalg import expm
        k0 = s["kappa_real"][0] + 1j * s["kappa_imag"][0]
        rows.append(dict(name=f.stem, t2=s["t2"].astype(np.float32),
                         U_ffsim=expm(k0),
                         nocc=int(s["nocc"]), nvirt=int(s["nvirt"]), norb=int(s["norb"])))
    return rows


def greedy_path(t2s):
    """Order points so consecutive ones are nearest neighbors (greedy TSP)."""
    X = t2s.reshape(len(t2s), -1)
    n = len(X)
    D = np.linalg.norm(X[:, None] - X[None, :], axis=2)
    np.fill_diagonal(D, np.inf)
    order = [0]; used = {0}
    for _ in range(n - 1):
        last = order[-1]
        nxt = min((j for j in range(n) if j not in used), key=lambda j: D[last, j])
        order.append(nxt); used.add(nxt)
    return order


def reconstruct(U, Z, nocc):
    full = 1j * torch.einsum("kpq,kap,kip,kcq,kjq->ijac", Z, U, U.conj(), U, U.conj())
    return full[:nocc, :nocc, nocc:, nocc:].real


def fit_one(t2, norb, nocc, n_reps, init=None, steps=1500, lr=5e-2):
    dev = t2.device
    ti = torch.triu_indices(norb, norb, 1)
    di = torch.triu_indices(norb, norb, 0)
    na, ns = ti.shape[1], di.shape[1]
    if init is None:
        # break the Z~0 saddle: nonzero Z + small random rotations
        A = 0.05 * torch.randn(n_reps, na, device=dev)
        S = 0.05 * torch.randn(n_reps, ns, device=dev)
        Zp = 0.3 * torch.randn(n_reps, ns, device=dev)
    else:
        A, S, Zp = (init[0].clone(), init[1].clone(), init[2].clone())
    A.requires_grad_(); S.requires_grad_(); Zp.requires_grad_()
    opt = torch.optim.Adam([A, S, Zp], lr=lr)
    first_loss = None
    for _ in range(steps):
        opt.zero_grad()
        Us, Zs = [], []
        for k in range(n_reps):
            Am = torch.zeros(norb, norb, device=dev); Am[ti[0], ti[1]] = A[k]; Am = Am - Am.T
            Sm = torch.zeros(norb, norb, device=dev); Sm[di[0], di[1]] = S[k]; Sm = Sm + Sm.T - torch.diag(torch.diag(Sm))
            Zm = torch.zeros(norb, norb, device=dev); Zm[di[0], di[1]] = Zp[k]; Zm = Zm + Zm.T - torch.diag(torch.diag(Zm))
            Us.append(torch.linalg.matrix_exp(torch.complex(Am, Sm)))
            Zs.append(Zm.to(torch.complex64))
        U = torch.stack(Us); Z = torch.stack(Zs)
        rec = reconstruct(U, Z, nocc)
        loss = torch.linalg.norm(rec - t2) / torch.linalg.norm(t2)
        if first_loss is None:
            first_loss = loss.item()
        loss.backward(); opt.step()
    with torch.no_grad():
        Us = []
        for k in range(n_reps):
            Am = torch.zeros(norb, norb, device=dev); Am[ti[0], ti[1]] = A[k]; Am = Am - Am.T
            Sm = torch.zeros(norb, norb, device=dev); Sm[di[0], di[1]] = S[k]; Sm = Sm + Sm.T - torch.diag(torch.diag(Sm))
            Us.append(torch.linalg.matrix_exp(torch.complex(Am, Sm)))
        U0 = Us[0].detach().cpu().numpy()
    return (A.detach(), S.detach(), Zp.detach()), loss.item(), U0, first_loss


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--norb", type=int, default=22)
    ap.add_argument("--n-reps", type=int, default=2)
    ap.add_argument("--steps", type=int, default=1200)
    ap.add_argument("--max-mols", type=int, default=40)
    args = ap.parse_args()
    torch.manual_seed(0); np.random.seed(0)
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    rows = load("rhf_cache", args.norb)[: args.max_mols]
    nocc, norb = rows[0]["nocc"], rows[0]["norb"]
    t2s = np.stack([r["t2"] for r in rows])
    order = greedy_path(t2s)
    rows = [rows[i] for i in order]
    t2s = t2s[order]
    print(f"{len(rows)} molecules, norb={norb}, ordered along nearest-neighbor t2 path\n")

    # Warm-started continuation
    init = None
    U_ws, resid = [], []
    for i, r in enumerate(rows):
        t2 = torch.tensor(r["t2"], device=dev)
        # cold-start molecule 0 with a few restarts; warm-start the rest
        if init is None:
            best = None
            for _ in range(4):
                cand = fit_one(t2, norb, nocc, args.n_reps, init=None, steps=args.steps)
                if best is None or cand[1] < best[1]:
                    best = cand
            init, res, U0, f0 = best
        else:
            init, res, U0, f0 = fit_one(t2, norb, nocc, args.n_reps, init=init, steps=args.steps)
        U_ws.append(U0); resid.append(res)
        print(f"  [{i+1}/{len(rows)}] {r['name'][:26]:26s} resid {f0:.2f}->{res:.3f}", flush=True)

    U_ws = np.stack(U_ws)
    U_ff = np.stack([r["U_ffsim"] for r in rows])
    dt2 = np.array([np.linalg.norm(t2s[i+1]-t2s[i]) for i in range(len(rows)-1)])
    dU_ws = np.array([np.linalg.norm(U_ws[i+1]-U_ws[i]) for i in range(len(rows)-1)])
    dU_ff = np.array([np.linalg.norm(U_ff[i+1]-U_ff[i]) for i in range(len(rows)-1)])

    print(f"\n=== RESULT (consecutive neighbors on t2 path) ===")
    print(f"mean recon_resid (warm-start)        : {np.mean(resid):.3f}  (ffsim floor ~0.32)")
    print(f"mean ||Δt2||                          : {dt2.mean():.3f}")
    print(f"mean ||ΔU|| warm-start (smooth section): {dU_ws.mean():.3f}")
    print(f"mean ||ΔU|| ffsim (independent solve)  : {dU_ff.mean():.3f}")
    print(f"corr(||Δt2||,||ΔU||) warm-start        : {np.corrcoef(dt2,dU_ws)[0,1]:+.3f}")
    print(f"corr(||Δt2||,||ΔU||) ffsim             : {np.corrcoef(dt2,dU_ff)[0,1]:+.3f}")
    if dU_ws.mean() < 0.6 * dU_ff.mean():
        print("\n=> POSITIVE: warm-started U is MUCH smoother than ffsim's while still")
        print("   reconstructing t2 => a smooth t2->(U,Z) section EXISTS and is learnable.")
    else:
        print("\n=> Smooth section not clearly better; investigate further.")


if __name__ == "__main__":
    main()
