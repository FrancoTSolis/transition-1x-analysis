#!/usr/bin/env python3
"""How good is a chemistry frame?  Residual of U = Phi B expm(A) with exact (VarPro) Z, for A = 0 and after
batched Adam on real generators A_k with growing support (same atom -> + bonded atoms -> unrestricted),
compared with the optimize=True labels (from the canonical exact-DF init) on the same molecules.

Usage: python3 -m pretrain.opt_true.frame_eval [--variant hyb_oao] [--split val] [--steps 300 300 400]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import varpro as V  # noqa: E402
from pretrain.opt_true.data import N29, load_split  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def real_sector_phases(n, nocc, device, dtype=torch.complex128):
    """(2, n) diagonal phases: rep 0 = (1_occ, e^{-i pi/4}), rep 1 = conj."""
    ph = torch.ones(n, dtype=dtype, device=device)
    ph[nocc:] = complex(math.cos(math.pi / 4), -math.sin(math.pi / 4))
    return torch.stack([ph, ph.conj()])


def u_from_real(B, A, phases):
    """B (b,2,n,n) real orthogonal, A (b,2,n,n) real (antisymmetrized here) -> complex U (b,2,n,n)."""
    K = A - A.transpose(-1, -2)
    O = B @ torch.linalg.matrix_exp(K)
    return phases[None, :, :, None] * O.to(phases.dtype)


def load_frames(names, device, dtype=torch.float64, shape="16_13"):
    f = np.load(ROOT / "rhf_frames" / f"{shape}.npz")
    pos = {str(k): i for i, k in enumerate(f["names"])}
    sel = [pos[k] for k in names]
    keys = [k for k in ("B", "slot_atom", "slot_type", "slot_partner", "slot_dir", "atom_Z", "coords", "bond", "natm")
            if k in f.files]
    out = {k: f[k][sel] for k in keys}
    t = {"B": torch.as_tensor(out["B"], device=device).to(dtype)}
    sa = torch.as_tensor(out["slot_atom"].astype(np.int64), device=device)
    bond = torch.as_tensor(out["bond"], device=device)
    b = torch.arange(len(sel), device=device)[:, None, None]
    same = (sa[:, :, None] == sa[:, None, :])
    bonded = bond[b, sa[:, :, None], sa[:, None, :]] | same
    if "slot_partner" in out:
        t["slot_partner"] = torch.as_tensor(out["slot_partner"].astype(np.int64), device=device)
        t["slot_dir"] = torch.as_tensor(out["slot_dir"], device=device)
    t["bond_full"] = bond
    t.update(slot_atom=sa, same=same, bonded=bonded, slot_type=torch.as_tensor(out["slot_type"].astype(np.int64), device=device),
             atom_Z=torch.as_tensor(out["atom_Z"].astype(np.int64), device=device),
             coords=torch.as_tensor(out["coords"], device=device).to(dtype), natm=out["natm"])
    return t


N_CHEM_NODE = 4 + 5 + 5 + 3 + 3 + 1
N_CHEM_PAIR = 1 + 1 + 8 + 3 + 2 + 3


def chem_features(fr):
    """Static chemistry features of the frame slots: node (b,n,21), pair (b,n,n,18), float32."""
    sa, st = fr["slot_atom"], fr["slot_type"]
    b = torch.arange(sa.shape[0], device=sa.device)[:, None]
    Zs = fr["atom_Z"][b, sa]                                                   # (b,n)
    el = torch.stack([Zs == z for z in (1, 6, 7, 8)], -1).float()
    t0 = torch.nn.functional.one_hot(st[:, 0], 10)[..., :5].float()
    t1 = torch.nn.functional.one_hot(st[:, 1], 10)[..., [0, 6, 7, 8, 9]].float()
    sdir = fr["slot_dir"].float()
    xyz = fr["coords"][b, sa].float()                                          # (b,n,3)
    deg = fr["bond_full"].sum(-1).float()[b, sa] / 4.0
    node = torch.cat([el, t0, t1, sdir, xyz / 3.0, deg[..., None]], -1)
    rel = xyz[:, None, :, :] - xyz[:, :, None, :]                               # r_q - r_p
    dist = rel.norm(dim=-1)
    cent = torch.linspace(0.0, 5.0, 8, device=sa.device)
    rbf = torch.exp(-((dist[..., None] - cent) / 0.7) ** 2)
    part = fr["slot_partner"]
    p_to_q = (part[:, :, None] == sa[:, None, :]).float()
    q_to_p = (part[:, None, :] == sa[:, :, None]).float()
    u = rel / dist.clamp_min(1e-6)[..., None]
    dd = (sdir[:, :, None, :] * sdir[:, None, :, :]).sum(-1)
    dpu = (sdir[:, :, None, :] * u).sum(-1)
    dqu = -(sdir[:, None, :, :] * u).sum(-1)
    pair = torch.cat([fr["same"].float()[..., None], fr["bonded"].float()[..., None], rbf, rel / 3.0,
                      p_to_q[..., None], q_to_p[..., None], dd[..., None], dpu[..., None], dqu[..., None]], -1)
    return node, pair


def variant_B(fr, variant):
    B = fr["B"]
    if variant == "hyb_oao":
        return B
    if variant == "oao_oao":
        return torch.stack([B[:, 1], B[:, 1]], 1)
    if variant == "hyb_hyb":
        return torch.stack([B[:, 0], B[:, 0]], 1)
    if variant == "oao_hyb":
        return torch.stack([B[:, 1], B[:, 0]], 1)
    raise ValueError(variant)


def varpro_eval(U, t2, mask, lam, zref):
    Z = V.solve_z(U, t2, mask, lam, zref)
    return Z, D.rel_residual(Z, U, t2), D.rel_residual_imag(Z, U, t2), D.normalized_objective(Z, U, t2, lam, zref)


def optimize(B, phases, t2, mask, lam, zref, supports, steps, lr=0.02, A0=None, log_every=0):
    A = (torch.zeros_like(B) if A0 is None else A0.clone()).requires_grad_(True)
    hist = []
    for sup, T in zip(supports, steps):
        opt = torch.optim.Adam([A], lr=lr)
        for it in range(T):
            U = u_from_real(B, A * sup, phases)
            with torch.no_grad():
                Z = V.solve_z(U, t2, mask, lam, zref)
            F = D.normalized_objective(Z, U, t2, lam, zref)
            opt.zero_grad(set_to_none=True)
            F.sum().backward()
            opt.step()
            if log_every and (it + 1) % log_every == 0:
                hist.append(float(F.median()))
        with torch.no_grad():
            A.mul_(sup)
        U = u_from_real(B, A.detach(), phases)
        Z, r, ri, F = varpro_eval(U, t2, mask, lam, zref)
        hist.append((float(r.median()), float(ri.median()), float(F.median())))
    return A.detach(), hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="hyb_oao")
    ap.add_argument("--split", default="val")
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--steps", type=int, nargs=3, default=[300, 300, 400])
    ap.add_argument("--lr", type=float, default=0.02)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--chunk", type=int, default=72)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    dev = "cuda"
    tr, va = load_split()
    names = va if args.split == "val" else tr
    if args.n:
        names = names[:args.n]
    d = N29(names=names, device=dev, dtype=torch.float64)
    fr = load_frames(names, dev)
    mask = D.square_mask(29, dev, torch.float64)
    ph = real_sector_phases(29, 16, dev)
    B = variant_B(fr, args.variant)
    t0 = time.time()
    res = {"variant": args.variant, "n": len(names)}
    rows = {k: [] for k in ("frame", "same_atom", "bonded", "full")}
    lab, canon = [], []
    for s in range(0, len(names), args.chunk):
        sl = slice(s, s + args.chunk)
        t2, zr = d.t2[sl], d.znorm_full[sl]
        U = u_from_real(B[sl], torch.zeros_like(B[sl]), ph)
        _, r, ri, F = varpro_eval(U, t2, mask, args.lam, zr)
        rows["frame"].append(torch.stack([r, ri, F], 1))
        Zl = d.Z_opt[sl]
        lab.append(torch.stack([D.rel_residual(Zl, d.U_opt[sl], t2), D.rel_residual_imag(Zl, d.U_opt[sl], t2),
                                D.normalized_objective(Zl, d.U_opt[sl], t2, args.lam, zr)], 1))
        _, rc, rci, Fc = varpro_eval(d.U0[sl], t2, mask, args.lam, zr)
        canon.append(torch.stack([rc, rci, Fc], 1))
        same = fr["same"][sl][:, None].to(B.dtype).expand(-1, 2, -1, -1)
        bonded = fr["bonded"][sl][:, None].to(B.dtype).expand(-1, 2, -1, -1)
        full = torch.ones_like(same)
        A = None
        for tag, sup, T in (("same_atom", same, args.steps[0]), ("bonded", bonded, args.steps[1]),
                            ("full", full, args.steps[2])):
            A, _ = optimize(B[sl], ph, t2, mask, args.lam, zr, [sup], [T], lr=args.lr, A0=A)
            U = u_from_real(B[sl], A, ph)
            _, r, ri, F = varpro_eval(U, t2, mask, args.lam, zr)
            rows[tag].append(torch.stack([r, ri, F], 1))
        print(f"  chunk {s}: {time.time()-t0:.0f}s", flush=True)
    med = lambda x: [round(float(v), 4) for v in torch.cat(x).median(0).values]  # noqa: E731
    res["label(optimize=True from canonical)"] = med(lab)
    res["canonical init + VarPro Z"] = med(canon)
    for k, v in rows.items():
        res[f"{args.variant}:{k}"] = med(v)
    print(json.dumps(res, indent=1))
    print("columns: [real residual, imag residual, normalized objective] (medians)")
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
