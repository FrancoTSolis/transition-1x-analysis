#!/usr/bin/env python3
"""Train the slot-frame learned optimizer (pretrain.opt_true.slotnet) on the n29 group, objective only.

Loss: mean over recycles t=1..T of the normalized ffsim objective F(U_t, Z*(U_t)) (Z* detached = Danskin
gradient), with weight 2 on the last recycle. T_train is sampled per batch from [1, --T] (so the model is
usable at any recycle count). Eval: val/train median residual (and objective) after t = 0, 1, 2, 4, 8, ... T_eval.

Usage: python3 -m pretrain.opt_true.train_slot --out runs_ot/slot_T1 --T 1 --epochs 60
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
from pretrain.opt_true.data import N29, load_split  # noqa: E402
from pretrain.opt_true.slotnet import SlotNet  # noqa: E402
from pretrain.opt_true import frame_eval as FE  # noqa: E402


def evaluate(model, d, idx, mask, lam, T_eval, bs=48, ctx=None):
    model.eval()
    res = {}
    for s in range(0, len(idx), bs):
        b = idx[s:s + bs]
        with torch.no_grad():
            traj = model(ctx["U0"][b], d.t2[b], mask, lam, d.znorm_full[b], T=T_eval,
                         chem=ctx["chem"](b), gen_mask=ctx["gm"](b))
        for t, (U, Z, obj) in enumerate(traj):
            r = D.rel_residual(Z, U, d.t2[b])
            ri = D.rel_residual_imag(Z, U, d.t2[b])
            res.setdefault(t, {"r": [], "i": [], "o": []})
            res[t]["r"].append(r); res[t]["i"].append(ri); res[t]["o"].append(obj)
    model.train()
    out = {}
    for t, v in res.items():
        out[t] = {"resid": float(torch.cat(v["r"]).median()), "imag": float(torch.cat(v["i"]).median()),
                  "obj": float(torch.cat(v["o"]).median())}
    return out


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--T", type=int, default=1, help="max recycles during training")
    ap.add_argument("--T-eval", type=int, default=0, help="recycles at eval (default: T, and 2T if T>1)")
    ap.add_argument("--fixed-T", action="store_true", help="always train with exactly T recycles")
    ap.add_argument("--d", type=int, default=128)
    ap.add_argument("--layers", type=int, default=4)
    ap.add_argument("--heads", type=int, default=8)
    ap.add_argument("--kscale", type=float, default=1.0)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--warmup", type=int, default=200)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--bs", type=int, default=16)
    ap.add_argument("--clip", type=float, default=1.0)
    ap.add_argument("--train-subset", type=int, default=0)
    ap.add_argument("--eval-every", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--frame", default="canonical", choices=["canonical", "chem"])
    ap.add_argument("--variant", default="hyb_oao", help="chem frame variant (frame_eval.variant_B)")
    ap.add_argument("--no-chem-feats", action="store_true")
    ap.add_argument("--gen-support", default="full", choices=["full", "bonded", "same"])
    ap.add_argument("--init-from", default=None)
    args = ap.parse_args()
    torch.manual_seed(args.seed)
    dev = "cuda"
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    (out / "args.json").write_text(json.dumps(vars(args), indent=1))
    logf = open(out / "log.jsonl", "a")

    def log(rec):
        rec["time"] = time.time(); logf.write(json.dumps(rec) + "\n"); logf.flush()
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

    d = N29(config=None, device=dev)
    mask = D.square_mask(29, dev)
    tr, va = load_split()
    itr, iva = d.indices(tr), d.indices(va)
    if args.train_subset:
        itr = itr[:args.train_subset]
    T_eval = args.T_eval or (args.T if args.T == 1 else 2 * args.T)
    chem_on = args.frame == "chem" and not args.no_chem_feats
    if args.frame == "chem":
        fr = FE.load_frames(d.names, dev)
        ph = FE.real_sector_phases(29, 16, dev, torch.complex64)
        Bf = FE.variant_B(fr, args.variant).float()
        U0all = ph[None, :, :, None] * Bf.to(torch.complex64)
        cnode, cpair = FE.chem_features(fr) if chem_on else (None, None)
        gsup = {"full": None, "bonded": fr["bonded"], "same": fr["same"]}[args.gen_support]
    else:
        U0all, cnode, cpair, gsup = d.U0, None, None, None
    ctx = {"U0": U0all,
           "chem": (lambda b: (cnode[b], cpair[b])) if chem_on else (lambda b: None),
           "gm": (lambda b: gsup[b]) if gsup is not None else (lambda b: None)}
    model = SlotNet(d=args.d, layers=args.layers, heads=args.heads, kscale=args.kscale,
                    max_recycles=max(16, T_eval + 1), real=args.frame == "chem",
                    n_chem_node=FE.N_CHEM_NODE if chem_on else 0, n_chem_pair=FE.N_CHEM_PAIR if chem_on else 0).to(dev)
    if args.init_from:
        sd = torch.load(args.init_from, map_location=dev, weights_only=False)
        print("init:", model.load_state_dict(sd["model"], strict=False))
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.0)
    spe = math.ceil(len(itr) / args.bs); total = args.epochs * spe
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / args.warmup) * 0.5 *
                                              (1 + math.cos(math.pi * min(1.0, s / max(1, total)))))
    log({"type": "start", "n_params": sum(p.numel() for p in model.parameters()), "n_train": len(itr), "steps": total})

    def do_eval(ep):
        rec = {"type": "eval", "epoch": ep}
        ev = evaluate(model, d, iva, mask, args.lam, T_eval, ctx=ctx)
        et = evaluate(model, d, itr[:141], mask, args.lam, T_eval, ctx=ctx)
        for t in sorted(ev):
            if t in (0, 1, 2, 4, 8, 16, 32, 64) or t == T_eval or t == args.T:
                rec[f"val_t{t}"] = ev[t]["resid"]; rec[f"val_imag_t{t}"] = ev[t]["imag"]; rec[f"val_obj_t{t}"] = ev[t]["obj"]
                rec[f"train_t{t}"] = et[t]["resid"]
        return rec, ev[args.T]["resid"]

    rec, best = do_eval(0); log(rec)
    g = torch.Generator().manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    for ep in range(1, args.epochs + 1):
        t0 = time.time(); tot = 0.0; cnt = 0
        perm = itr[torch.randperm(len(itr), generator=g).to(dev)]
        for s in range(0, len(perm), args.bs):
            b = perm[s:s + args.bs]
            T = args.T if (args.fixed_T or args.T == 1) else int(rng.integers(1, args.T + 1))
            traj = model(ctx["U0"][b], d.t2[b], mask, args.lam, d.znorm_full[b], T=T, chem=ctx["chem"](b),
                         gen_mask=ctx["gm"](b))
            objs = [D.normalized_objective(Z.detach(), U, d.t2[b], args.lam, d.znorm_full[b]) for (U, Z, _) in traj[1:]]
            w = torch.ones(len(objs), device=dev); w[-1] = 2.0
            loss = sum(wi * o.mean() for wi, o in zip(w, objs)) / w.sum()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), args.clip)
            opt.step(); sched.step()
            tot += float(objs[-1].mean()) * len(b); cnt += len(b)
        rec = {"type": "train", "epoch": ep, "last_obj": tot / cnt, "time": time.time() - t0, "lr": sched.get_last_lr()[0]}
        if ep % args.eval_every == 0 or ep == args.epochs:
            r, v = do_eval(ep); rec.update(r)
            if v < best:
                best = v
                torch.save({"model": model.state_dict(), "args": vars(args), "epoch": ep}, out / "best.pt")
        log(rec)
    torch.save({"model": model.state_dict(), "args": vars(args), "epoch": args.epochs}, out / "last.pt")
    log({"type": "done", "best_val": best})


if __name__ == "__main__":
    main()
