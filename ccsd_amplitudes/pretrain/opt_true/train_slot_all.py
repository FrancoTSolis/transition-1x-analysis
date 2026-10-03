#!/usr/bin/env python3
"""Multi-size training of the chemistry-frame slot model (pretrain.opt_true.slotnet, real sector) on every
molecule shape, objective only (exact VarPro Z, Danskin gradient).  Batches never mix shapes; batch size per
shape follows a budget on B * n^3 (edge-transformer cost).  Held out: the n29 val split (141) and the 30 small
molecules (norb <= 16) used for exact energies.  Eval: n29 val / n29 train subset / small molecules at
t = 0, 1, T (and 2T if T > 1).

Usage: python3 -m pretrain.opt_true.train_slot_all --out runs_ot/slotall_T1 --T 1 --epochs 6
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
from pretrain.opt_true import frame_eval as FE  # noqa: E402
from pretrain.opt_true.data import load_split  # noqa: E402
from pretrain.opt_true.slotnet import SlotNet  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FRAME_KEYS = ("B", "slot_atom", "slot_type", "slot_partner", "slot_dir", "atom_Z", "coords", "bond")


class Bucket:
    """All molecules of one shape (CPU tensors); .batch(ids, dev) returns GPU tensors + chemistry features."""

    def __init__(self, no, nv, names, variant="hyb_oao"):
        self.no, self.nv, self.n = no, nv, no + nv
        f = np.load(ROOT / "rhf_frames" / f"{no}_{nv}.npz")
        pos = {str(k): i for i, k in enumerate(f["names"])}
        ci = np.load(ROOT / "rhf_inits_canonical" / f"{no}_{nv}.npz")
        cpos = {str(k): i for i, k in enumerate(ci["names"])}
        names = [k for k in names if k in pos and k in cpos]
        self.names, self.N = names, len(names)
        sel = [pos[k] for k in names]
        self.fr = {k: torch.as_tensor(f[k][sel]) for k in FRAME_KEYS}
        self.variant = variant
        self.znorm = torch.as_tensor(ci["znorm_full"][[cpos[k] for k in names]]).float()
        self.t2 = torch.as_tensor(np.stack([np.load(ROOT / "rhf_dataset" / f"{k}.npz")["t2"] for k in names])).float()

    def batch(self, ids, dev):
        ids = torch.as_tensor(ids)
        fr = {k: v[ids].to(dev) for k, v in self.fr.items()}
        b = torch.arange(len(ids), device=dev)[:, None, None]
        sa = fr["slot_atom"].long()
        fr["slot_atom"], fr["slot_type"], fr["slot_partner"] = sa, fr["slot_type"].long(), fr["slot_partner"].long()
        fr["atom_Z"] = fr["atom_Z"].long()
        same = sa[:, :, None] == sa[:, None, :]
        fr["same"], fr["bond_full"] = same, fr["bond"]
        fr["bonded"] = fr["bond"][b, sa[:, :, None], sa[:, None, :]] | same
        fr["coords"] = fr["coords"].float()
        B = FE.variant_B(fr, self.variant).float()
        ph = FE.real_sector_phases(self.n, self.no, dev, torch.complex64)
        U0 = ph[None, :, :, None] * B.to(torch.complex64)
        node, pair = FE.chem_features(fr)
        return {"U0": U0, "t2": self.t2[ids].to(dev), "zref": self.znorm[ids].to(dev), "chem": (node, pair),
                "mask": D.square_mask(self.n, dev)}


def evaluate(model, buckets, lam, T_eval, dev, bs=32):
    """Median residual at each recycle over the union of the given buckets."""
    model.eval()
    res = {}
    for bk in buckets:
        for s in range(0, bk.N, bs):
            x = bk.batch(list(range(s, min(s + bs, bk.N))), dev)
            with torch.no_grad():
                traj = model(x["U0"], x["t2"], x["mask"], lam, x["zref"], T=T_eval, chem=x["chem"])
            for t, (U, Z, obj) in enumerate(traj):
                res.setdefault(t, []).append(D.rel_residual(Z, U, x["t2"]))
    model.train()
    return {t: float(torch.cat(v).median()) for t, v in res.items()}


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--T", type=int, default=1)
    ap.add_argument("--fixed-T", action="store_true")
    ap.add_argument("--d", type=int, default=128)
    ap.add_argument("--layers", type=int, default=4)
    ap.add_argument("--heads", type=int, default=8)
    ap.add_argument("--kscale", type=float, default=1.0)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--warmup", type=int, default=500)
    ap.add_argument("--epochs", type=int, default=6)
    ap.add_argument("--budget", type=float, default=16 * 29 ** 3, help="sum over batch of n^3")
    ap.add_argument("--max-bs", type=int, default=64)
    ap.add_argument("--max-norb", type=int, default=44)
    ap.add_argument("--clip", type=float, default=1.0)
    ap.add_argument("--eval-every-steps", type=int, default=1500)
    ap.add_argument("--init-from", default=None)
    ap.add_argument("--variant", default="hyb_oao")
    ap.add_argument("--teacher", default=None, help="slot checkpoint whose T-recycle output the student imitates")
    ap.add_argument("--teacher-T", type=int, default=4)
    ap.add_argument("--w-distill", type=float, default=1.0, help="weight of ||U_student - U_teacher||^2 / (2n)")
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
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

    tr29, va29 = load_split()
    small = [ln.strip() for ln in open(ROOT / "gauge_study" / "names_small_norb_le16.txt") if ln.strip()]
    held = set(va29) | set(small)
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    by_shape = {}
    for k, (n, no, nv) in idx.items():
        if n <= args.max_norb and k not in held:
            by_shape.setdefault((no, nv), []).append(k)
    t0 = time.time()
    buckets = {s: Bucket(*s, sorted(v), args.variant) for s, v in by_shape.items()}
    buckets = {s: b for s, b in buckets.items() if b.N}
    keys = list(buckets)
    n_train = sum(b.N for b in buckets.values())
    bss = {s: int(max(1, min(args.max_bs, args.budget // (s[0] + s[1]) ** 3))) for s in keys}
    steps_per_epoch = int(sum(math.ceil(buckets[s].N / bss[s]) for s in keys))
    total = args.epochs * steps_per_epoch
    ev29 = Bucket(16, 13, va29, args.variant)
    ev29tr = Bucket(16, 13, tr29[:141], args.variant)
    small_shapes = {}
    for k in small:
        n, no, nv = idx[k]
        small_shapes.setdefault((no, nv), []).append(k)
    ev_small = [Bucket(no, nv, v, args.variant) for (no, nv), v in small_shapes.items()]
    T_eval = args.T if args.T == 1 else 2 * args.T
    from pretrain.opt_true.frame_eval import N_CHEM_NODE, N_CHEM_PAIR
    model = SlotNet(d=args.d, layers=args.layers, heads=args.heads, kscale=args.kscale, max_recycles=max(16, T_eval + 1),
                    real=True, n_chem_node=N_CHEM_NODE, n_chem_pair=N_CHEM_PAIR).to(dev)
    if args.init_from:
        sd = torch.load(args.init_from, map_location=dev, weights_only=False)
        print("init:", model.load_state_dict(sd["model"], strict=False))
    teacher = None
    if args.teacher:
        from pretrain.opt_true.eval_slot import load_slot
        teacher = load_slot(args.teacher, dev)[0]
        for p_ in teacher.parameters():
            p_.requires_grad_(False)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.0)
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1.0, (s + 1) / args.warmup) * 0.5 *
                                              (1 + math.cos(math.pi * min(1.0, s / max(1, total)))))
    log({"type": "start", "n_train": n_train, "n_buckets": len(keys), "steps_per_epoch": steps_per_epoch,
         "total_steps": total, "n_params": sum(p.numel() for p in model.parameters()), "load_s": time.time() - t0})

    def do_eval(step):
        rec = {"type": "eval", "step": step}
        for tag, bks in (("val29", [ev29]), ("train29", [ev29tr]), ("small", ev_small)):
            r = evaluate(model, bks, args.lam, T_eval, dev)
            for t, v in r.items():
                if t in (0, 1, args.T, T_eval):
                    rec[f"{tag}_t{t}"] = v
        return rec

    best, step = float("inf"), 0
    log(do_eval(0))
    for ep in range(1, args.epochs + 1):
        plan = []
        for s in keys:
            perm = rng.permutation(buckets[s].N)
            plan += [(s, perm[i:i + bss[s]]) for i in range(0, len(perm), bss[s])]
        rng.shuffle(plan)
        t1, tot, cnt = time.time(), 0.0, 0
        for s, ids in plan:
            x = buckets[s].batch(ids, dev)
            T = args.T if (args.fixed_T or args.T == 1) else int(rng.integers(1, args.T + 1))
            traj = model(x["U0"], x["t2"], x["mask"], args.lam, x["zref"], T=T, chem=x["chem"])
            objs = [D.normalized_objective(Z.detach(), U, x["t2"], args.lam, x["zref"]) for (U, Z, _) in traj[1:]]
            w = [1.0] * (len(objs) - 1) + [2.0]
            loss = sum(wi * o.mean() for wi, o in zip(w, objs)) / sum(w)
            if teacher is not None:
                with torch.no_grad():
                    Ut = teacher(x["U0"], x["t2"], x["mask"], args.lam, x["zref"], T=args.teacher_T, chem=x["chem"])[-1][0]
                Us = traj[-1][0]
                ld = (Us - Ut).abs().pow(2).flatten(1).sum(1) / (2 * Us.shape[-1])
                loss = loss + args.w_distill * ld.mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), args.clip)
            opt.step()
            sched.step()
            step += 1
            tot += float(objs[-1].mean()) * len(ids)
            cnt += len(ids)
            if step % args.eval_every_steps == 0:
                rec = do_eval(step)
                rec.update(epoch=ep, train_obj_running=tot / max(cnt, 1))
                log(rec)
                v = rec[f"val29_t{args.T}"]
                if v < best:
                    best = v
                    torch.save({"model": model.state_dict(), "args": vars(args), "step": step}, out / "best.pt")
        log({"type": "epoch", "epoch": ep, "train_obj": tot / max(cnt, 1), "time": time.time() - t1, "step": step})
        torch.save({"model": model.state_dict(), "args": vars(args), "step": step}, out / "last.pt")
    rec = do_eval(step)
    log(rec)
    if rec[f"val29_t{args.T}"] < best:
        torch.save({"model": model.state_dict(), "args": vars(args), "step": step}, out / "best.pt")
    log({"type": "done", "best_val29": min(best, rec[f"val29_t{args.T}"])})


if __name__ == "__main__":
    main()
