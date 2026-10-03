#!/usr/bin/env python3
"""Multi-size optimize=True pretraining: every molecule shape in rhf_dataset/, objective-only (no labels).

Same model / parameterization / unrolled steps as pretrain.opt_true.train (canonical frame, right-multiplied
generator, pair Z, optional --unroll T), but data are bucketed by exact shape (nocc, nvirt) so a batch never
needs padding. Batch size per bucket follows a token budget (the positional slice encoder scales as
n^2 * nocc * nvirt). Held out from training: the n29 val split (141) and the 30 small molecules (norb <= 16)
used for exact-energy evaluation. Reported each eval: n29 val, n29 train subset and small-molecule residuals.

Usage: python3 -m pretrain.opt_true.train_all --out runs_ot/all_unroll8 --unroll 8 --epochs 20
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
from pretrain.model import ModelConfig  # noqa: E402
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import train as T  # noqa: E402
from pretrain.opt_true.data import load_split  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


class Bucket:
    """All molecules of one shape, CPU-resident; mimics the N29 API used by train.evaluate."""

    def __init__(self, no, nv, names, init_npz, dev):
        self.no, self.nv, self.n = no, nv, no + nv
        pos = {str(k): i for i, k in enumerate(init_npz["names"])}
        sel = [pos[k] for k in names]
        self.names = names
        self.N = len(names)
        self.t2 = torch.as_tensor(np.stack([np.load(ROOT / "rhf_dataset" / f"{k}.npz")["t2"] for k in names])).float()
        self.U0 = torch.complex(torch.as_tensor(init_npz["U_re"][sel]), torch.as_tensor(init_npz["U_im"][sel]))
        self.Z0_full = torch.as_tensor(init_npz["Z"][sel]).float()
        self.znorm_full = torch.as_tensor(init_npz["znorm_full"][sel]).float()
        self.mask = D.square_mask(self.n)
        self.dev = dev

    def gpu(self):
        """A view of this bucket with tensors on the GPU (for evaluate)."""
        g = Bucket.__new__(Bucket)
        g.__dict__.update(self.__dict__)
        g.t2, g.U0 = self.t2.to(self.dev), self.U0.to(self.dev)
        g.Z0_full, g.znorm_full = self.Z0_full.to(self.dev), self.znorm_full.to(self.dev)
        g.mask = self.mask.to(self.dev)
        return g

    def make_batch(self, idx):
        B = len(idx)
        dev = self.t2.device
        return {"t2": self.t2[idx], "noccs": torch.full((B,), self.no, device=dev, dtype=torch.long),
                "nvirts": torch.full((B,), self.nv, device=dev, dtype=torch.long),
                "norbs": torch.full((B,), self.n, device=dev, dtype=torch.long),
                "nocc": [self.no] * B, "nvirt": [self.nv] * B, "n_reps": [2] * B,
                "max_nocc": self.no, "max_nvirt": self.nv}


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--unroll", type=int, default=8)
    ap.add_argument("--unroll-eta", type=float, default=0.05)
    ap.add_argument("--unroll-deep", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--kscale", type=float, default=3.0)
    ap.add_argument("--zscale", type=float, default=1.0)
    ap.add_argument("--embed-dim", type=int, default=192)
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--heads", type=int, default=8)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--warmup", type=int, default=500)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--token-budget", type=float, default=12 * 29 ** 2 * 16 * 13, help="sum_batch n^2*no*nv")
    ap.add_argument("--max-norb", type=int, default=44)
    ap.add_argument("--clip", type=float, default=1.0)
    ap.add_argument("--eval-every-steps", type=int, default=2000)
    ap.add_argument("--init-from", default=None)
    ap.add_argument("--t2-scale", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    dev = "cuda"
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "args.json").write_text(json.dumps(vars(args), indent=1))
    logf = open(out / "log.jsonl", "a")

    def log(rec):
        rec["time"] = time.time()
        logf.write(json.dumps(rec) + "\n")
        logf.flush()
        print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

    tr29, va29 = load_split()
    small = [ln.strip() for ln in open(ROOT / "gauge_study" / "names_small_norb_le16.txt") if ln.strip()]
    held = set(va29) | set(small)
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    buckets, n_train = {}, 0
    for f in sorted((ROOT / "rhf_inits_canonical").glob("*.npz")):
        no, nv = map(int, f.stem.split("_"))
        if no + nv > args.max_norb:
            continue
        z = np.load(f)
        names = [str(k) for k in z["names"] if str(k) not in held and str(k) in idx]
        if names:
            buckets[(no, nv)] = Bucket(no, nv, names, z, dev)
            n_train += len(names)
    keys = list(buckets)
    sizes = np.array([buckets[k].N for k in keys], float)
    cost = np.array([(k[0] + k[1]) ** 2 * k[0] * k[1] for k in keys], float)
    bss = np.maximum(1, np.floor(args.token_budget / cost)).astype(int)
    steps_per_epoch = int(np.sum(np.ceil(sizes / bss)))
    total = args.epochs * steps_per_epoch
    # eval sets
    z29 = np.load(ROOT / "rhf_inits_canonical" / "16_13.npz")
    ev29 = Bucket(16, 13, va29, z29, dev).gpu()
    ev29tr = Bucket(16, 13, tr29[:141], z29, dev).gpu()
    small_by_shape = {}
    for k in small:
        n, no, nv = idx[k]
        small_by_shape.setdefault((no, nv), []).append(k)
    ev_small = [Bucket(no, nv, nm, np.load(ROOT / "rhf_inits_canonical" / f"{no}_{nv}.npz"), dev).gpu()
                for (no, nv), nm in small_by_shape.items()]

    mcfg = ModelConfig(embed_dim=args.embed_dim, num_layers=args.layers, num_heads=args.heads, n_reps=2, dropout=0.0,
                       attention_dropout=0.0, predict_residual=True, residual_kappa_scale=args.kscale,
                       residual_z_scale=args.zscale, residual_zero_init=True, residual_hyps=1)
    model = T.OTModel(mcfg, D.square_mask(29), "pair", input_scale=args.t2_scale).to(dev)
    unroller = T.Unroller(args.unroll, args.unroll_eta).to(dev) if args.unroll else None
    if args.init_from:
        sd = torch.load(args.init_from, map_location=dev, weights_only=False)
        print("init:", model.load_state_dict(sd["model"], strict=False))
        if unroller is not None and sd.get("unroller"):
            unroller.load_state_dict(sd["unroller"])
    params = list(model.parameters()) + (list(unroller.parameters()) if unroller else [])
    opt = torch.optim.AdamW(params, lr=args.lr, weight_decay=0.0)
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / args.warmup) * 0.5 *
                                              (1 + math.cos(math.pi * min(1.0, s / max(1, total)))))
    T.PARAM_MODE = "right"
    log({"type": "start", "n_train": n_train, "n_buckets": len(keys), "steps_per_epoch": steps_per_epoch,
         "total_steps": total, "n_params": sum(p.numel() for p in model.parameters())})

    def do_eval(step):
        rec = {"type": "eval", "step": step}
        for tag, bk in (("val29", ev29), ("train29", ev29tr)):
            r = T.evaluate(model, bk, torch.arange(bk.N, device=dev), bk.U0, bk.Z0_full * bk.mask, bk.mask,
                           args.lam, None, unroller=unroller)
            rec[f"{tag}_resid"] = r["resid_median"]
            rec[f"{tag}_imag"] = r["imag_median"]
            rec[f"{tag}_obj"] = r["obj_median"]
        rs = []
        for bk in ev_small:
            r = T.evaluate(model, bk, torch.arange(bk.N, device=dev), bk.U0, bk.Z0_full * bk.mask, bk.mask,
                           args.lam, None, unroller=unroller)
            rs.append((r["resid_median"], bk.N))
        rec["small_resid"] = float(np.average([a for a, _ in rs], weights=[b for _, b in rs]))
        return rec

    best = float("inf")
    step = 0
    log(do_eval(0))
    for ep in range(1, args.epochs + 1):
        plan = []
        for bi, k in enumerate(keys):
            perm = rng.permutation(buckets[k].N)
            for s in range(0, len(perm), bss[bi]):
                plan.append((k, perm[s:s + bss[bi]]))
        rng.shuffle(plan)
        t0, tot, cnt = time.time(), 0.0, 0
        for k, ids in plan:
            bk = buckets[k]
            ids_t = torch.as_tensor(ids)
            batch = {kk: (v.to(dev, non_blocking=True) if torch.is_tensor(v) else v)
                     for kk, v in bk.make_batch(ids_t).items()}      # CPU bucket -> GPU batch
            t2 = batch["t2"]
            U0 = bk.U0[ids_t].to(dev)
            mask = bk.mask.to(dev)
            Z0 = bk.Z0_full[ids_t].to(dev) * mask
            zr = bk.znorm_full[ids_t].to(dev)
            outp = model(batch)
            h = outp["residual_hyps"][0]
            p = [h["dkappa_re"], h["dkappa_im"], h["dz"]]
            deep = None
            if unroller is not None:
                p, traj = unroller(p, U0, Z0, mask, t2, args.lam, zr, track=True)
                deep = torch.stack(traj).mean(0) if traj else None
            U, Z = T.assemble(p, U0, Z0, mask)
            lm = T.per_mol_objective(U, Z, t2, args.lam, zr)
            if deep is not None and args.unroll_deep:
                lm = lm + args.unroll_deep * deep
            loss = lm.mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, args.clip)
            opt.step()
            sched.step()
            step += 1
            tot += float(loss) * len(ids)
            cnt += len(ids)
            if step % args.eval_every_steps == 0:
                rec = do_eval(step)
                rec["epoch"] = ep
                rec["train_loss_running"] = tot / max(cnt, 1)
                log(rec)
                if rec["val29_resid"] < best:
                    best = rec["val29_resid"]
                    torch.save({"model": model.state_dict(), "unroller": unroller.state_dict() if unroller else None,
                                "args": {**vars(args), "frame": "canonical", "param": "right", "zhead": "pair", "hyps": 1,
                                         "hyp_kick": 0.0, "t2_scale": args.t2_scale}, "step": step}, out / "best.pt")
        log({"type": "epoch", "epoch": ep, "loss": tot / max(cnt, 1), "time": time.time() - t0, "step": step})
        torch.save({"model": model.state_dict(), "unroller": unroller.state_dict() if unroller else None,
                    "args": {**vars(args), "frame": "canonical", "param": "right", "zhead": "pair", "hyps": 1,
                             "hyp_kick": 0.0, "t2_scale": args.t2_scale}, "step": step}, out / "last.pt")
    log({"type": "done", "best_val29_resid": best})


if __name__ == "__main__":
    main()
