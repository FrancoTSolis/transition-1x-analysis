#!/usr/bin/env python3
"""Experiment 14: are optimize=True labels smooth/learnable once aligned over the FULL symmetry group?

exp7 only removed per-column phases. The compressed-DF objective also has the discrete symmetries
swap (rep order), conj (U -> conj U, Z -> -Z), rev (index reversal, square mask) and gamma (sign
flip of virtual rows); the exact-DF init is a conj-swap-symmetric saddle, so the optimizer breaks
the symmetry by round-off and each label is a random orbit representative.

Measures, on the 1410 n29 molecules (labels: square lambda=0.005):
  A. closest-decile / all-pairs ratio of U distances: raw, phase-only, phase+subgroups, full group
  B. consistent representatives: iteratively align every label to a template (medoid), then
     ridge R^2 (t2-PCA256 features) of aligned U / kappa / Z vs the phase-only versions (exp7 recipe)
  C. orbit-distance 1-NN transfer: does the label of the t2-nearest train molecule predict the test
     label better than a random one, under the full-group distance?

Run (train venv, 1 GPU):  python3 -m gauge_study.exp14_orbit_aligned_labels --config square_reg0.005
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
from pretrain.opt_true import symmetry as G  # noqa: E402
from pretrain.opt_true.data import N29, load_split  # noqa: E402
from gauge_study.common import nearest_decile_ratio  # noqa: E402


def ridge_r2(Xtr, Ytr, Xte, Yte, alphas=(0.1, 1.0, 10.0, 100.0, 1000.0)):
    from sklearn.linear_model import Ridge
    from sklearn.metrics import r2_score
    n_val = max(1, int(0.15 * len(Xtr)))
    best, keep = -np.inf, 1.0
    for a in alphas:
        r = Ridge(a).fit(Xtr[:-n_val], Ytr[:-n_val])
        v = r2_score(Ytr[-n_val:].ravel(), r.predict(Xtr[-n_val:]).ravel())
        if v > best:
            best, keep = v, a
    r = Ridge(keep).fit(Xtr, Ytr)
    P = r.predict(Xte)
    mu = Ytr.mean(0, keepdims=True)
    pooled = r2_score(Yte.ravel(), P.ravel())
    spec = 1.0 - ((Yte - P) ** 2).sum() / max(((Yte - mu) ** 2).sum(), 1e-30)
    return pooled, spec, P


def herm_feats(U):
    """Phase-free per-column projector features: |U|^2 entries and U diag(w) U^dag (real/imag) with fixed weights."""
    n = U.shape[-1]
    w = torch.linspace(-1, 1, n, device=U.device, dtype=U.real.dtype).to(U.dtype)
    H = torch.einsum("brip,p,brjp->brij", U, w, U.conj())
    iu = torch.triu_indices(n, n)
    return torch.cat([H.real[..., iu[0], iu[1]].flatten(1), H.imag[..., iu[0], iu[1]].flatten(1)], 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="square_reg0.005")
    ap.add_argument("--out", default="gauge_study/exp14_results.json")
    args = ap.parse_args()
    t0 = time.time()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    d = N29(config=args.config, device=dev, dtype=torch.float64)
    N = d.N
    U, Z = d.U_opt, d.Z_opt
    t2f = d.t2.flatten(1)
    sq = (t2f ** 2).sum(1)
    Dt2 = (sq[:, None] + sq[None] - 2 * t2f @ t2f.T).clamp_min(0).sqrt()
    iu = torch.triu_indices(N, N, 1, device=dev)
    d_in = Dt2[iu[0], iu[1]].cpu().numpy()
    res = {"config": args.config, "N": N}
    print(f"{N} molecules ({time.time()-t0:.0f}s)", flush=True)

    # ---------------- A. smoothness under increasing alignment ----------------
    print("\n--- A. closest-decile / all-pairs ratio of label U distances (input t2 ratio %.3f) ---"
          % nearest_decile_ratio(d_in, d_in), flush=True)
    raw = (U[:, None] - U[None]).abs().pow(2).flatten(2).sum(2) if N <= 0 else None
    groups = {
        "phase only": [(0, 0, 0, 0)],
        "phase+swap": G.elements(True, False, False, False),
        "phase+conj": G.elements(False, True, False, False),
        "phase+swap+conj": G.elements(True, True, False, False),
        "phase+swap+conj+rev": G.elements(True, True, True, False),
        "phase+full (swap,conj,rev,gamma)": G.elements(True, True, True, True),
    }
    res["ratios"] = {}
    for name, grp in groups.items():
        Dm = G.pairwise_orbit_sqdist(U, grp).clamp_min(0).sqrt()
        v = Dm[iu[0], iu[1]].cpu().numpy()
        r = nearest_decile_ratio(d_in, v)
        c = float(np.corrcoef(d_in, v)[0, 1])
        res["ratios"][name] = {"ratio": float(r), "corr": c, "median_dist": float(np.median(v))}
        print(f"  {name:36s} ratio={r:.3f} corr={c:+.3f} median pair dist={np.median(v):.3f}", flush=True)
        if name.startswith("phase+full"):
            Dfull = Dm

    # how often does the best pairwise group element differ from identity among nearest pairs?
    # ---------------- B. consistent representatives + ridge ----------------
    print("\n--- B. align all labels to a template (iterated medoid), ridge R^2 ---", flush=True)
    tr, va = load_split()
    itr, iva = d.indices(tr), d.indices(va)
    med = int(torch.argmin(Dfull[itr][:, itr].sum(1)))
    ref_U, ref_Z = U[itr][med:med + 1], Z[itr][med:med + 1]
    for it in range(4):
        Ua, Za, gidx, dist = G.orbit_align(ref_U.expand_as(U), ref_Z.expand_as(Z), U, Z)
        # new template = (phase-insensitive) mean of aligned train labels -> nearest label
        mean_U = Ua[itr].mean(0, keepdim=True)
        cand = G.phase_sqdist(mean_U.expand_as(Ua[itr]), Ua[itr])
        k = int(torch.argmin(cand))
        ref_U, ref_Z = Ua[itr][k:k + 1], Za[itr][k:k + 1]
        counts = torch.bincount(gidx, minlength=16).tolist()
        print(f"  iter {it}: mean dist to template {dist.sqrt().mean():.3f}; group-element usage {counts}", flush=True)
    res["group_usage"] = counts
    from sklearn.decomposition import PCA
    X = d.t2.flatten(1).cpu().numpy().astype(np.float32)
    pca = PCA(n_components=256, svd_solver="randomized", random_state=0).fit(X[itr.cpu().numpy()])
    Xtr, Xte = pca.transform(X[itr.cpu().numpy()]), pca.transform(X[iva.cpu().numpy()])

    def feats_U(Ux):
        return torch.cat([Ux.real.flatten(1), Ux.imag.flatten(1)], 1).cpu().numpy()

    # phase-only alignment to the same template (exp7-like but via template phases)
    Up, Zp, _, _ = G.orbit_align(ref_U.expand_as(U), ref_Z.expand_as(Z), U, Z, group=[(0, 0, 0, 0)])
    res["ridge"] = {}
    for name, Ux, Zx in [("phase-only aligned", Up, Zp), ("full-group aligned", Ua, Za)]:
        Yu = feats_U(Ux)
        pooled, spec, P = ridge_r2(Xtr, Yu[itr.cpu()], Xte, Yu[iva.cpu()])
        Yz = Zx.flatten(1).cpu().numpy()
        keep = np.abs(Yz).sum(0) > 0
        zp, zs, _ = ridge_r2(Xtr, Yz[itr.cpu()][:, keep], Xte, Yz[iva.cpu()][:, keep])
        # decode U prediction -> nearest unitary (polar) -> phase dist to truth, relative
        Pc = torch.as_tensor(P[:, :Yu.shape[1] // 2] + 1j * P[:, Yu.shape[1] // 2:], device=dev).reshape(-1, 2, 29, 29)
        Uu, _, Vh = torch.linalg.svd(Pc)
        Upol = Uu @ Vh
        du = (G.phase_sqdist(Ux[iva], Upol).sqrt() / np.sqrt(58)).median().item()
        du_mean = (G.phase_sqdist(Ux[iva], Ux[itr].mean(0, keepdim=True).expand_as(Ux[iva])).sqrt() / np.sqrt(58)).median().item()
        res["ridge"][name] = {"U_pooled": pooled, "U_specific": float(spec), "Z_pooled": zp, "Z_specific": float(zs),
                              "U_rel_dist_pred": du, "U_rel_dist_trainmean": du_mean}
        print(f"  {name:20s} U: R2 pooled {pooled:+.3f} specific {spec:+.3f} | Z: pooled {zp:+.3f} specific {zs:+.3f} | "
              f"val U rel dist pred {du:.3f} vs train-mean {du_mean:.3f}", flush=True)

    # ---------------- C. 1-NN transfer under full-group distance ----------------
    Dtt = Dt2[iva][:, itr]
    nn = torch.argmin(Dtt, 1)
    rnd = torch.randint(0, len(itr), (len(iva),), device=dev)
    Dv = Dfull[iva][:, itr]
    r1 = Dv[torch.arange(len(iva)), nn].mean().item()
    rr = Dv[torch.arange(len(iva)), rnd].mean().item()
    res["nn_transfer_full"] = r1 / rr
    print(f"\n--- C. 1-NN label transfer, full-group distance: 1-NN {r1:.3f} vs random {rr:.3f} (ratio {r1/rr:.3f}) ---")
    json.dump(res, open(args.out, "w"), indent=1)
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
