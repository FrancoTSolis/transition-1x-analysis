#!/usr/bin/env python3
"""Attribute the gain of an unrolled model: network vs learned optimizer steps vs plain optimization.

On the n29 val split, for a checkpoint with an unroller (T learned steps):
  init              canonical exact-DF init
  net               network output only (no steps)
  zero+learnedT     the SAME learned T steps started from the init (network output replaced by 0)
  net+learnedT      the full model
  init+adamT        plain Adam (lr 0.01) for T steps from the init      (same step budget)
  init+adam{N}      plain Adam for N steps from the init (reference optimum, ~optimize=True)
  net+learnedT+adam{N}  the full model followed by N Adam steps
Reports median real residual, imaginary residual and normalized ffsim objective.

Usage: python3 -m pretrain.opt_true.decompose --ckpt runs_ot/unroll32/best.pt [--adam-ref 500]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import train as T  # noqa: E402
from pretrain.opt_true.data import N29, load_split  # noqa: E402
from pretrain.opt_true.eval_energy import load_model  # noqa: E402


def stats(U, Z, t2, lam, zr):
    return {"resid": float(D.rel_residual(Z, U, t2).median()), "imag": float(D.rel_residual_imag(Z, U, t2).median()),
            "obj": float(D.normalized_objective(Z, U, t2, lam, zr).median())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--adam-ref", type=int, default=500)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    dev = "cuda"
    model, unroller, a = load_model(args.ckpt, dev)
    Tn = unroller.T if unroller is not None else 0
    _, va = load_split()
    d = N29(names=va, config=None, device=dev)
    m = D.square_mask(29, dev)
    Z0 = d.Z0_full * m
    res = {}
    U0, t2, zr = d.U0, d.t2, d.znorm_full
    res["init"] = stats(U0, Z0, t2, args.lam, zr)
    outs = []
    with torch.no_grad():
        for s in range(0, d.N, 24):
            b = torch.arange(s, min(s + 24, d.N), device=dev)
            h = model(d.make_batch(b))["residual_hyps"][0]
            outs.append([h["dkappa_re"], h["dkappa_im"], h["dz"]])
    p_net = [torch.cat([o[i] for o in outs]) for i in range(3)]
    U, Z = T.assemble(p_net, U0, Z0, m)
    res["net"] = stats(U, Z, t2, args.lam, zr)

    def unroll_from(p):
        outs2 = []
        for s in range(0, d.N, 24):
            sl = slice(s, s + 24)
            with torch.enable_grad():
                pp = [x[sl].detach().clone().requires_grad_(True) for x in p]
                pT, _ = unroller(pp, U0[sl], Z0[sl], m, t2[sl], args.lam, zr[sl])
            outs2.append([x.detach() for x in pT])
        return [torch.cat([o[i] for o in outs2]) for i in range(3)]

    if unroller is not None:
        unroller.eval()
        p_zero_T = unroll_from([torch.zeros_like(x) for x in p_net])
        res[f"zero+learned{Tn}"] = stats(*T.assemble(p_zero_T, U0, Z0, m), t2, args.lam, zr)
        p_full = unroll_from(p_net)
        res[f"net+learned{Tn}"] = stats(*T.assemble(p_full, U0, Z0, m), t2, args.lam, zr)
        Uf, Zf, _ = D.refine(t2, U0, Z0, m, lam=args.lam, znorm_ref=zr, steps=args.adam_ref, lr=0.01,
                             A0=p_full[0], S0=p_full[1], D0=p_full[2])
        res[f"net+learned{Tn}+adam{args.adam_ref}"] = stats(Uf, Zf, t2, args.lam, zr)
        Ua, Za, _ = D.refine(t2, U0, Z0, m, lam=args.lam, znorm_ref=zr, steps=Tn, lr=0.01)
        res[f"init+adam{Tn}"] = stats(Ua, Za, t2, args.lam, zr)
        for lr in (0.03, 0.1):
            Ua, Za, _ = D.refine(t2, U0, Z0, m, lam=args.lam, znorm_ref=zr, steps=Tn, lr=lr)
            res[f"init+adam{Tn}(lr{lr})"] = stats(Ua, Za, t2, args.lam, zr)
    Ur, Zr, _ = D.refine(t2, U0, Z0, m, lam=args.lam, znorm_ref=zr, steps=args.adam_ref, lr=0.01)
    res[f"init+adam{args.adam_ref}"] = stats(Ur, Zr, t2, args.lam, zr)
    print(f"checkpoint {args.ckpt} (T={Tn}, epoch {torch.load(args.ckpt, map_location='cpu', weights_only=False).get('epoch')}) on {d.N} n29 val molecules:")
    for k, v in res.items():
        print(f"  {k:30s} resid {v['resid']:.3f}  imag {v['imag']:.3f}  obj {v['obj']:.3f}")
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
