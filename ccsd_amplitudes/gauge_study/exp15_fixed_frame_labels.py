#!/usr/bin/env python3
"""Experiment 15: optimize=True labels from a MOLECULE-INDEPENDENT starting frame.

The canonical exact-DF init U0(t2) is essentially random across molecules (its columns live in
near-degenerate eigen-pairs), and it is a conj-swap-symmetric saddle that the optimizer leaves in
a round-off-chosen direction. Both make the labels a random function of t2.

Here every molecule starts from the SAME frame: U_k = U_ref,k (fixed seeded random unitaries,
different for the two reps, so the conjugate-pair symmetry is broken identically everywhere) and
Z = 0, and ffsim's complex objective (lambda=0.005, square mask) is minimized with batched Adam on
all 1410 molecules at once. If the optimizer's solution depends smoothly on t2 when the start is
shared, the labels become consistent:
  - quality: median relative residual vs the canonical labels (0.342 on val, 0.328 all)
  - smoothness: closest-decile ratio of raw (no alignment needed: shared frame) U distances, t2hat
  - learnability: ridge R^2 (t2-PCA256) of U entries / Z, molecule-specific
Variants: --frame random (Haar, seeded) | randsmall (expm of small random generator) ; --zinit 0|lam

Run (train venv, 1 GPU):  python3 -m gauge_study.exp15_fixed_frame_labels --frame random
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gauge_study.common import nearest_decile_ratio  # noqa: E402
from gauge_study.exp14_orbit_aligned_labels import ridge_r2  # noqa: E402
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true.data import N29, load_split  # noqa: E402


def haar_unitary(n, gen, dtype=torch.complex128):
    X = torch.complex(torch.randn(n, n, generator=gen, dtype=torch.float64),
                      torch.randn(n, n, generator=gen, dtype=torch.float64))
    Q, R = torch.linalg.qr(X)
    ph = torch.diagonal(R) / torch.diagonal(R).abs()
    return (Q * ph[None, :]).to(dtype)


def complex_objective_refine(t2, Ub, Zb, mask, lam, zref, steps, lr, log_every=500):
    """Adam on ffsim's complex objective; U = Ub expm(K), Z = (Zb + D) * mask."""
    B, R, n, _ = Ub.shape
    rdt = t2.dtype
    A = torch.zeros(B, R, n, n, device=t2.device, dtype=rdt, requires_grad=True)
    S = torch.zeros_like(A, requires_grad=True)
    Dp = torch.zeros_like(A, requires_grad=True)
    opt = torch.optim.Adam([A, S, Dp], lr=lr, betas=(0.9, 0.99))
    t2n2 = t2.flatten(1).pow(2).sum(1)
    hist = []
    for it in range(steps + 1):
        U = D.unitary_from(Ub, A, S)
        Z = (Zb + D.sym(Dp)) * mask
        rec = D.reconstruct_complex(Z, U, t2.shape[1])
        r = rec - t2
        loss = 0.5 * r.abs().pow(2).flatten(1).sum(1) + lam * (Z.flatten(1).pow(2).sum(1) - zref).abs()
        if it % log_every == 0 or it == steps:
            with torch.no_grad():
                rr = r.real.flatten(1).norm(dim=1) / t2n2.sqrt()
                ri = r.imag.flatten(1).norm(dim=1) / t2n2.sqrt()
                hist.append((it, float(rr.median()), float(ri.median())))
                print(f"    step {it:5d}: resid_real median {rr.median():.4f}  resid_imag {ri.median():.4f}", flush=True)
        if it == steps:
            break
        opt.zero_grad(set_to_none=True)
        (loss / (0.5 * t2n2)).sum().backward()
        opt.step()
    with torch.no_grad():
        U = D.unitary_from(Ub, A, S)
        Z = (Zb + D.sym(Dp)) * mask
    return U.detach(), Z.detach(), hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", default="random", choices=["random", "randsmall"])
    ap.add_argument("--scale", type=float, default=0.5, help="[randsmall] generator scale")
    ap.add_argument("--zinit", default="0", choices=["0", "small"])
    ap.add_argument("--steps", type=int, default=3000)
    ap.add_argument("--lr", type=float, default=0.01)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=470)
    ap.add_argument("--save", default=None, help="write labels npz to this path")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    t0 = time.time()
    dev = "cuda"
    dt = torch.float32
    d = N29(config="square_reg0.005", device=dev, dtype=dt)
    N, n = d.N, 29
    m = D.square_mask(n, dev, dt)
    gen = torch.Generator().manual_seed(args.seed)
    if args.frame == "random":
        Uref = torch.stack([haar_unitary(n, gen), haar_unitary(n, gen)])
    else:
        Ks = []
        for _ in range(2):
            X = torch.complex(torch.randn(n, n, generator=gen, dtype=torch.float64),
                              torch.randn(n, n, generator=gen, dtype=torch.float64))
            Ks.append(args.scale * (X - X.conj().T) / (2 * np.sqrt(n)))
        Uref = torch.linalg.matrix_exp(torch.stack(Ks))
    Uref = Uref.to(dev).to(torch.complex64)
    Ub = Uref[None].expand(N, 2, n, n).contiguous()
    if args.zinit == "0":
        Zb = torch.zeros(N, 2, n, n, device=dev, dtype=dt)
    else:
        g2 = torch.Generator().manual_seed(args.seed + 1)
        Zb = (0.05 * torch.randn(2, n, n, generator=g2)).to(dev).to(dt)
        Zb = (D.sym(Zb) * m)[None].expand(N, 2, n, n).contiguous()
    print(f"frame={args.frame} zinit={args.zinit}: optimizing {N} molecules ({time.time()-t0:.0f}s)", flush=True)
    Us, Zs, hists = [], [], []
    for s in range(0, N, args.chunk):
        sl = slice(s, s + args.chunk)
        print(f"  chunk {s}:{s+args.chunk}", flush=True)
        U, Z, h = complex_objective_refine(d.t2[sl], Ub[sl], Zb[sl], m, args.lam, d.znorm_full[sl], args.steps, args.lr)
        Us.append(U)
        Zs.append(Z)
        hists.append(h)
    U = torch.cat(Us).to(torch.complex128)
    Z = torch.cat(Zs).to(torch.float64)
    t2 = d.t2.to(torch.float64)
    if args.save:   # save before the analysis so a later failure cannot lose the labels
        np.savez(args.save, names=np.array(d.names), U=U.cpu().numpy(), Z=Z.cpu().numpy(), Uref=Uref.cpu().numpy(),
                 resid=D.rel_residual(Z, U, t2).cpu().numpy())
    rr = D.rel_residual(Z, U, t2)
    ri = (D.reconstruct_complex(Z, U, 16).imag.flatten(1).norm(dim=1) / t2.flatten(1).norm(dim=1))
    tr, va = load_split()
    itr, iva = d.indices(tr), d.indices(va)
    res = {"args": vars(args), "resid_real_median": float(rr.median()), "resid_real_val_median": float(rr[iva].median()),
           "resid_imag_median": float(ri.median()), "label_canonical_median": float(d.resid_opt.median()),
           "label_canonical_val_median": float(d.resid_opt[iva].median()), "znorm_median": float(Z.flatten(1).norm(dim=1).median())}
    print(f"\nquality: fixed-frame resid median {rr.median():.4f} (val {rr[iva].median():.4f}) imag {ri.median():.4f} | "
          f"canonical labels {d.resid_opt.median():.4f} (val {d.resid_opt[iva].median():.4f})", flush=True)

    # smoothness (shared frame -> raw distances are meaningful; also phase-aligned)
    t2f = t2.flatten(1)
    sq = (t2f ** 2).sum(1)
    Dt2 = (sq[:, None] + sq[None] - 2 * t2f @ t2f.T).clamp_min(0).sqrt()
    iu = torch.triu_indices(N, N, 1, device=dev)
    d_in = Dt2[iu[0], iu[1]].cpu().numpy()
    Uf = U.flatten(1)
    Dr = ((Uf.abs() ** 2).sum(1)[:, None] + (Uf.abs() ** 2).sum(1)[None] - 2 * (Uf.conj() @ Uf.T).real).clamp_min(0).sqrt()
    from pretrain.opt_true.symmetry import pairwise_orbit_sqdist
    Dp_ = pairwise_orbit_sqdist(U, [(0, 0, 0, 0)]).clamp_min(0).sqrt()
    Zf = Z.flatten(1)
    Dz = torch.cdist(Zf, Zf)
    th = D.reconstruct(Z, U, 16).flatten(1)
    Dth = torch.cdist(th, th)
    res["ratios"] = {}
    for name, M in [("U raw (shared frame)", Dr), ("U phase-aligned", Dp_), ("Z", Dz), ("t2hat", Dth)]:
        v = M[iu[0], iu[1]].cpu().numpy()
        res["ratios"][name] = float(nearest_decile_ratio(d_in, v))
        print(f"  smoothness {name:22s} ratio={res['ratios'][name]:.3f} corr={np.corrcoef(d_in, v)[0,1]:+.3f}", flush=True)

    # learnability
    from sklearn.decomposition import PCA
    X = t2f.cpu().numpy().astype(np.float32)
    itr_np, iva_np = itr.cpu().numpy(), iva.cpu().numpy()
    pca = PCA(n_components=256, svd_solver="randomized", random_state=0).fit(X[itr_np])
    Xtr, Xte = pca.transform(X[itr_np]), pca.transform(X[iva_np])
    Yu = torch.cat([U.real.flatten(1), U.imag.flatten(1)], 1).cpu().numpy()
    pu, su, P = ridge_r2(Xtr, Yu[itr_np], Xte, Yu[iva_np])
    Yz = Z.flatten(1).cpu().numpy()
    keep = np.abs(Yz).sum(0) > 0
    pz, sz, Pz = ridge_r2(Xtr, Yz[itr_np][:, keep], Xte, Yz[iva_np][:, keep])
    # decode ridge prediction -> nearest unitary; evaluate the residual of the PREDICTED (U, Z)
    h = Yu.shape[1] // 2
    Pc = torch.as_tensor(P[:, :h] + 1j * P[:, h:], device=dev).reshape(-1, 2, n, n)
    W, _, Vh = torch.linalg.svd(Pc)
    Upred = W @ Vh
    Zpred = torch.zeros(len(iva_np), 2 * n * n, device=dev, dtype=torch.float64)
    Zpred[:, torch.as_tensor(keep, device=dev)] = torch.as_tensor(Pz, device=dev).double()
    Zpred = D.sym(Zpred.reshape(-1, 2, n, n)) * m.double()
    r_pred = D.rel_residual(Zpred, Upred, t2[iva])
    r_mean = D.rel_residual(Z[itr].mean(0, keepdim=True).expand(len(iva_np), 2, n, n),
                            Upred[:1].expand(len(iva_np), 2, n, n) * 0 + U[itr][:1].expand(len(iva_np), 2, n, n), t2[iva])
    res["ridge"] = {"U_pooled": pu, "U_specific": float(su), "Z_pooled": pz, "Z_specific": float(sz),
                    "val_resid_of_ridge_prediction": float(r_pred.median())}
    print(f"  ridge: U R2 pooled {pu:+.3f} specific {su:+.3f} | Z pooled {pz:+.3f} specific {sz:+.3f} | "
          f"val residual of ridge-predicted (U,Z): {r_pred.median():.3f}", flush=True)
    out = args.out or f"gauge_study/exp15_{args.frame}_z{args.zinit}_s{args.seed}.json"
    json.dump(res, open(out, "w"), indent=1)
    print(f"total {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
