#!/usr/bin/env python3
"""Evaluate slot-model checkpoints (train_slot / train_slot_all) on the n29 splits.

Splits: val141 (the fixed val split), val89 (leak-free: no train molecule with the same HF energy; the 52 removed
are all reactants with bit- or near-identical train twins) and train141 (first 141 train molecules).
Reports median real residual, imaginary residual and normalized ffsim objective after t = 0, 1, 2, 4, 8 recycles,
next to the optimize=True labels (from the canonical init) on the same molecules.

Usage: python3 -m pretrain.opt_true.eval_slot --ckpt runs_ot/slotchem2_T1/best.pt [--ckpt ...] [--T-max 8]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import frame_eval as FE  # noqa: E402
from pretrain.opt_true.data import N29, load_split  # noqa: E402
from pretrain.opt_true.slotnet import SlotNet  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def load_slot(ckpt, dev):
    sd = torch.load(ckpt, map_location=dev, weights_only=False)
    a = sd["args"]
    chem = a.get("frame", "chem") == "chem" and not a.get("no_chem_feats", False)
    real = a.get("frame", "chem") == "chem"
    m = SlotNet(d=a["d"], layers=a["layers"], heads=a["heads"], kscale=a.get("kscale", 1.0),
                max_recycles=sd["model"]["t_emb.weight"].shape[0], real=real,
                n_chem_node=FE.N_CHEM_NODE if chem else 0, n_chem_pair=FE.N_CHEM_PAIR if chem else 0).to(dev)
    m.load_state_dict(sd["model"])
    m.eval()
    return m, a, chem, real


def context(d, a, chem, dev):
    if a.get("frame", "chem") == "chem":
        fr = FE.load_frames(d.names, dev)
        ph = FE.real_sector_phases(29, 16, dev, torch.complex64)
        U0 = ph[None, :, :, None] * FE.variant_B(fr, a.get("variant", "hyb_oao")).float().to(torch.complex64)
        cn, cp = FE.chem_features(fr) if chem else (None, None)
        gs = {"full": None, "bonded": fr["bonded"], "same": fr["same"]}[a.get("gen_support", "full")]
    else:
        U0, cn, cp, gs = d.U0, None, None, None
    return U0, cn, cp, gs


@torch.no_grad()
def run(model, d, idx, U0, cn, cp, gs, lam, T, bs=48):
    out = {}
    for s in range(0, len(idx), bs):
        b = idx[s:s + bs]
        traj = model(U0[b], d.t2[b], D.square_mask(29, d.t2.device), lam, d.znorm_full[b], T=T,
                     chem=(cn[b], cp[b]) if cn is not None else None, gen_mask=gs[b] if gs is not None else None)
        for t, (U, Z, obj) in enumerate(traj):
            o = out.setdefault(t, {"r": [], "i": [], "o": []})
            o["r"].append(D.rel_residual(Z, U, d.t2[b]))
            o["i"].append(D.rel_residual_imag(Z, U, d.t2[b]))
            o["o"].append(obj)
    return {t: {k: torch.cat(v[k]) for k in v} for t, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", action="append", required=True)
    ap.add_argument("--T-max", type=int, default=8)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    dev = "cuda"
    tr, va = load_split()
    clean = [ln.strip() for ln in open(ROOT / "pretrain" / "opt_true" / "split_n29_val_clean.txt") if ln.strip()]
    names = va + tr[:141]
    d = N29(names=names, device=dev)
    sets = {"val141": d.indices(va), "val89": d.indices(clean), "train141": d.indices(tr[:141])}
    lab = {k: (D.rel_residual(d.Z_opt[i], d.U_opt[i], d.t2[i]).median().item(),
               D.normalized_objective(d.Z_opt[i], d.U_opt[i], d.t2[i], args.lam, d.znorm_full[i]).median().item())
           for k, i in sets.items()}
    res = {"labels(optimize=True)": {k: {"resid": round(v[0], 4), "obj": round(v[1], 4)} for k, v in lab.items()}}
    for ck in args.ckpt:
        model, a, chem, real = load_slot(ck, dev)
        U0, cn, cp, gs = context(d, a, chem, dev)
        allidx = torch.arange(d.N, device=dev)
        r = run(model, d, allidx, U0, cn, cp, gs, args.lam, args.T_max)
        rr = {}
        for t in sorted(r):
            if t not in (0, 1, 2, 4, 8, 16, 32) and t != args.T_max:
                continue
            rr[f"t{t}"] = {k: {"resid": round(r[t]["r"][i].median().item(), 4), "imag": round(r[t]["i"][i].median().item(), 4),
                               "obj": round(r[t]["o"][i].median().item(), 4)} for k, i in sets.items()}
        res[ck] = {"train_T": a.get("T"), "frame": a.get("frame", "chem"), **rr}
    for k, v in res.items():
        print(k)
        if k.startswith("labels"):
            print("   " + "  ".join(f"{s}: {x['resid']:.3f}" for s, x in v.items()))
            continue
        for t, x in v.items():
            if t[0] == "t" and t[1:].isdigit():
                print(f"   {t:4s} " + "  ".join(f"{s}: {y['resid']:.3f} (obj {y['obj']:.3f})" for s, y in x.items()))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        json.dump(res, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
