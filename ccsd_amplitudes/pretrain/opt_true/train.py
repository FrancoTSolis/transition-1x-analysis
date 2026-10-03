#!/usr/bin/env python3
"""Train the Edge Transformer to output optimize=True-quality LUCJ parameters for the n29 group.

    tokens = PretrainingModel(t2)            (pair tokens -> residual heads, K hypotheses)
    --param right: U_k = U_base,k @ expm(dK_k)   dK acts on the COLUMNS (slots) of U_base
    --param left : U_k = expm(dK_k) @ U_base,k   dK acts on the ROWS (MO index) -- matches the
                   MO-indexed pair tokens: dK[p,q] really mixes orbitals p and q
    dK anti-Hermitian (tanh-bounded entries, --kscale)
    Z_k    = (Z_base,k + dZ_k) * mask        square connectivity; dZ from the pair heads (--zhead pair,
             slot-indexed outputs read from MO-indexed tokens) or from a global pooled head (--zhead global)

Frames (--frame):
    canonical : U_base = canonical exact-DF init U0(t2), Z_base = its masked Z (old compressed_recon)
    fixed     : U_base = a molecule-independent unitary per rep (from --frame-file, e.g. the U_ref of
                gauge_study/exp15), Z_base = 0 (or --zbase file)
Losses (weights may be combined):
    --w-obj   ffsim compressed-DF objective, COMPLEX reconstruction (|rec - t2|^2 / ||t2||^2 + 2 lam|.|/||t2||^2)
    --w-sup   supervised to labels: phase-invariant U loss mean_p(1-|<u_p,v_p>|^2) + --w-supz * ||Z-Zl||^2/||Zl||^2
              labels: --labels canonical:<config> (rhf_targets_compressed) or npz:<path> (exp15 format)
              --orbit-sup  takes the min of the supervised loss over the 16-element discrete symmetry group
    with K hypotheses the per-molecule loss is 0.95*min_k + 0.05*mean_k; inference picks the hypothesis with
    the lowest objective (free to evaluate).
Overfitting / capacity tests: --train-subset N trains on the first N train molecules and also reports the
train metrics on exactly those molecules each epoch.

Logs JSON lines to <out>/log.jsonl; saves <out>/best.pt (by val residual) and <out>/last.pt.
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
from pretrain.model import ModelConfig, PretrainingModel  # noqa: E402
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import symmetry as G  # noqa: E402
from pretrain.opt_true.data import N29, load_split  # noqa: E402


class GlobalZHead(torch.nn.Module):
    """dZ (masked entries, both reps) from pooled final pair tokens; zero-init output."""

    def __init__(self, d: int, mask: torch.Tensor, n_reps: int = 2, hidden: int = 256):
        super().__init__()
        iu = torch.nonzero(torch.triu(mask) > 0)
        self.register_buffer("ii", iu[:, 0])
        self.register_buffer("jj", iu[:, 1])
        self.n_reps, self.n, self.m = n_reps, mask.shape[0], iu.shape[0]
        self.mlp = torch.nn.Sequential(torch.nn.Linear(2 * d, hidden), torch.nn.GELU(),
                                       torch.nn.Linear(hidden, hidden), torch.nn.GELU(),
                                       torch.nn.Linear(hidden, n_reps * self.m))
        torch.nn.init.zeros_(self.mlp[-1].weight)
        torch.nn.init.zeros_(self.mlp[-1].bias)

    def forward(self, x):
        B = x.shape[0]
        pooled = torch.cat([x.mean((1, 2)), x.diagonal(dim1=1, dim2=2).mean(-1)], -1)
        v = self.mlp(pooled).view(B, self.n_reps, self.m)
        Dm = x.new_zeros(B, self.n_reps, self.n, self.n)
        Dm[:, :, self.ii, self.jj] = v
        return Dm + Dm.transpose(-1, -2) - torch.diag_embed(torch.diagonal(Dm, dim1=-2, dim2=-1))


class OTModel(torch.nn.Module):
    """PretrainingModel backbone + residual heads (+ optional global Z heads, one per hypothesis)."""

    def __init__(self, mcfg, mask, zhead: str = "pair", input_scale: float = 1.0):
        super().__init__()
        self.net = PretrainingModel(mcfg)
        self.zhead = zhead
        # t2 entries are ~2e-3 RMS, ~10x below the positional embeddings they are mixed with in the slice
        # encoder; scaling the NETWORK INPUT (not the objective) makes the molecule-specific signal visible.
        self.input_scale = input_scale
        self._x = None
        self.net.final_norm.register_forward_hook(lambda mod, inp, out: setattr(self, "_x", out))
        if zhead == "global":
            self.zheads = torch.nn.ModuleList([GlobalZHead(mcfg.embed_dim, mask) for _ in range(max(1, mcfg.residual_hyps))])

    def forward(self, batch):
        if self.input_scale != 1.0:
            batch = {**batch, "t2": batch["t2"] * self.input_scale}
        out = self.net(batch)
        hyps = out.get("residual_hyps") or [out]
        if self.zhead == "global":
            hyps = [{**h, "dz": self.zheads[k](self._x)} for k, h in enumerate(hyps)]
        return {"residual_hyps": hyps}


PARAM_MODE = "right"


class Unroller(torch.nn.Module):
    """T learned, differentiable normalized-gradient steps on the per-molecule objective, applied to the
    network's raw outputs p0 = (dK_re, dK_im, dZ): p_{t+1} = p_t - eta_{t,g} * g / ||g||_molecule.
    Trained end to end (create_graph=True) -- the model is 'network + T learned optimizer steps'."""

    def __init__(self, T: int, eta0: float):
        super().__init__()
        self.T = T
        self.log_eta = torch.nn.Parameter(torch.full((max(T, 1), 3), math.log(eta0)))

    def forward(self, p, U_base, Z_base, mask, t2, lam, zref, track=False):
        traj = []
        for t in range(self.T):
            U, Z = assemble(p, U_base, Z_base, mask)
            L = per_mol_objective(U, Z, t2, lam, zref)
            if track:
                traj.append(L)
            g = torch.autograd.grad(L.sum(), p, create_graph=self.training)
            p = [x - torch.exp(self.log_eta[t, i]) * gi / (gi.flatten(1).norm(dim=1).view(-1, 1, 1, 1) + 1e-8)
                 for i, (x, gi) in enumerate(zip(p, g))]
        return p, traj


def assemble(p, U_base, Z_base, mask):
    kr, ki, dz = p
    E = torch.linalg.matrix_exp(torch.complex(kr, ki).to(U_base.dtype))
    U = U_base @ E if PARAM_MODE == "right" else E @ U_base
    return U, (Z_base + dz) * mask


def build_params(out, U_base, Z_base, mask, hyp_list=None):
    """Return list over hypotheses of (U, Z) for the batch."""
    hyps = hyp_list if hyp_list is not None else (out.get("residual_hyps") or [out])
    res = []
    for h in hyps:
        K = torch.complex(h["dkappa_re"], h["dkappa_im"]).to(U_base.dtype)   # antisym + i*sym by construction
        E = torch.linalg.matrix_exp(K)
        U = U_base @ E if PARAM_MODE == "right" else E @ U_base
        Z = (Z_base + h["dz"]) * mask
        res.append((U, Z))
    return res


def per_mol_objective(U, Z, t2, lam, zref):
    return D.normalized_objective(Z, U, t2, lam, zref, complex_obj=True)


def sup_loss(U, Z, Ul, Zl, w_z, orbit=False):
    if not orbit:
        lu = D.phase_invariant_u_loss(U, Ul)
        lz = (Z - Zl).flatten(1).pow(2).sum(1) / Zl.flatten(1).pow(2).sum(1).clamp_min(1e-12)
        return lu + w_z * lz
    best = None
    for g in G.elements():
        Ug, Zg = G.apply(g, Ul, Zl)
        lu = D.phase_invariant_u_loss(U, Ug)
        lz = (Z - Zg).flatten(1).pow(2).sum(1) / Zg.flatten(1).pow(2).sum(1).clamp_min(1e-12)
        v = lu + w_z * lz
        best = v if best is None else torch.minimum(best, v)
    return best


@torch.no_grad()
def evaluate(model, d, idx, U_base_all, Z_base_all, mask, lam, labels, bs=24, unroller=None):
    model.eval()
    if unroller is not None:
        unroller.eval()
    rr, ri, ob, du, best_k = [], [], [], [], []
    for s in range(0, len(idx), bs):
        b = idx[s:s + bs]
        out = model(d.make_batch(b))
        if unroller is not None and unroller.T > 0:
            h = out["residual_hyps"][0]
            with torch.enable_grad():
                p0 = [h["dkappa_re"].detach().requires_grad_(True), h["dkappa_im"].detach().requires_grad_(True),
                      h["dz"].detach().requires_grad_(True)]
                pT, _ = unroller(p0, U_base_all[b], Z_base_all[b], mask, d.t2[b], lam, d.znorm_full[b])
            cands = [tuple(x.detach() for x in assemble(pT, U_base_all[b], Z_base_all[b], mask))]
        else:
            cands = build_params(out, U_base_all[b], Z_base_all[b], mask)
        objs = torch.stack([per_mol_objective(U, Z, d.t2[b], lam, d.znorm_full[b]) for U, Z in cands])  # (K,B)
        kbest = objs.argmin(0)
        U = torch.stack([cands[k][0][i] for i, k in enumerate(kbest.tolist())])
        Z = torch.stack([cands[k][1][i] for i, k in enumerate(kbest.tolist())])
        rr.append(D.rel_residual(Z, U, d.t2[b]))
        ri.append(D.rel_residual_imag(Z, U, d.t2[b]))
        ob.append(objs.min(0).values)
        best_k.append(kbest)
        if labels is not None:
            du.append(D.phase_invariant_u_loss(U, labels[0][b]))
    cat = lambda x: torch.cat(x)  # noqa: E731
    r = {"resid_median": float(cat(rr).median()), "resid_mean": float(cat(rr).mean()),
         "imag_median": float(cat(ri).median()), "obj_median": float(cat(ob).median())}
    if du:
        r["u_loss_to_label_median"] = float(cat(du).median())
    if len(best_k) and cat(best_k).max() > 0:
        r["hyp_usage"] = torch.bincount(cat(best_k)).tolist()
    model.train()
    if unroller is not None:
        unroller.train()
    return r


@torch.no_grad()
def distill_targets(model, d, idx, U_base_all, Z_base_all, mask, lam, steps, lr, prox, kscale, bs=96):
    """Refine the network's own predictions per molecule (batched Adam, complex objective) -> targets
    (dK_re, dK_im, dZ) in the network's own coordinates (right-multiplied generator)."""
    model.eval()
    A_, S_, D_, res = [], [], [], []
    for s in range(0, len(idx), bs):
        b = idx[s:s + bs]
        h = model(d.make_batch(b))["residual_hyps"][0]
        with torch.enable_grad():
            U, Z, _, (A, S, Dz) = D.refine(d.t2[b], U_base_all[b], Z_base_all[b], mask, lam=lam,
                                           znorm_ref=d.znorm_full[b], steps=steps, lr=lr, prox_mu=prox,
                                           A0=h["dkappa_re"], S0=h["dkappa_im"], D0=h["dz"], return_params=True)
        A = 0.5 * (A - A.transpose(-1, -2))
        S = 0.5 * (S + S.transpose(-1, -2))
        Dz = 0.5 * (Dz + Dz.transpose(-1, -2)) * mask
        lim = 0.98 * kscale
        A_.append(A.clamp(-lim, lim)); S_.append(S.clamp(-lim, lim)); D_.append(Dz)
        res.append(D.rel_residual(Z, U, d.t2[b]))
    model.train()
    return torch.cat(A_), torch.cat(S_), torch.cat(D_), torch.cat(res)


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--frame", default="fixed", choices=["fixed", "canonical"])
    ap.add_argument("--param", default="right", choices=["right", "left"])
    ap.add_argument("--zhead", default="pair", choices=["pair", "global"])
    ap.add_argument("--frame-file", default="gauge_study/exp15_labels_random.npz",
                    help="npz with Uref (2,n,n) for --frame fixed")
    ap.add_argument("--labels", default=None, help="canonical:<config> | npz:<path>")
    ap.add_argument("--w-obj", type=float, default=1.0)
    ap.add_argument("--w-sup", type=float, default=0.0)
    ap.add_argument("--w-supz", type=float, default=1.0)
    ap.add_argument("--orbit-sup", action="store_true")
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--hyps", type=int, default=1)
    ap.add_argument("--hyp-kick", type=float, default=0.0)
    ap.add_argument("--kscale", type=float, default=3.0, help="tanh bound on dK entries")
    ap.add_argument("--zscale", type=float, default=1.0, help="tanh bound on dZ entries")
    ap.add_argument("--no-zero-init", action="store_true")
    ap.add_argument("--embed-dim", type=int, default=192)
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--heads", type=int, default=8)
    ap.add_argument("--dropout", type=float, default=0.0)
    ap.add_argument("--wd", type=float, default=0.0)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--warmup", type=int, default=200)
    ap.add_argument("--schedule", default="cosine", choices=["cosine", "constant"])
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--bs", type=int, default=12)
    ap.add_argument("--clip", type=float, default=1.0)
    ap.add_argument("--train-subset", type=int, default=0)
    ap.add_argument("--eval-every", type=int, default=1)
    ap.add_argument("--init-from", default=None)
    ap.add_argument("--t2-scale", type=float, default=1.0, help="multiply the network's t2 INPUT (objective unchanged)")
    ap.add_argument("--unroll", type=int, default=0, help="T learned differentiable optimizer steps inside the model")
    ap.add_argument("--unroll-eta", type=float, default=0.05)
    ap.add_argument("--unroll-deep", type=float, default=0.1, help="weight of intermediate-step objectives")
    ap.add_argument("--distill-every", type=int, default=0, help="expert iteration: refresh refined targets every E epochs")
    ap.add_argument("--distill-start", type=int, default=1)
    ap.add_argument("--distill-steps", type=int, default=200)
    ap.add_argument("--distill-lr", type=float, default=0.01)
    ap.add_argument("--distill-prox", type=float, default=0.0)
    ap.add_argument("--w-distill", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    global PARAM_MODE
    PARAM_MODE = args.param
    dev = "cuda"
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "args.json").write_text(json.dumps(vars(args), indent=1))
    logf = open(out_dir / "log.jsonl", "a")

    def log(rec):
        rec["time"] = time.time()
        logf.write(json.dumps(rec) + "\n")
        logf.flush()
        print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

    cfg_lab = None
    if args.labels and args.labels.startswith("canonical:"):
        cfg_lab = args.labels.split(":", 1)[1]
    d = N29(config=cfg_lab or "square_reg0.005", device=dev, dtype=torch.float32)
    n = 29
    mask = D.square_mask(n, dev, torch.float32)
    tr, va = load_split()
    itr, iva = d.indices(tr), d.indices(va)
    if args.train_subset:
        itr = itr[: args.train_subset]

    if args.frame == "canonical":
        U_base_all = d.U0
        Z_base_all = d.Z0_full * mask
    else:
        f = np.load(args.frame_file)
        Uref = torch.as_tensor(f["Uref"], device=dev).to(torch.complex64)
        U_base_all = Uref[None].expand(d.N, 2, n, n).contiguous()
        Z_base_all = torch.zeros(d.N, 2, n, n, device=dev)

    labels = None
    if args.labels:
        if cfg_lab:
            labels = (d.U_opt, d.Z_opt)
        else:
            f = np.load(args.labels.split(":", 1)[1], allow_pickle=True)
            pos = {str(nm): i for i, nm in enumerate(f["names"])}
            order = [pos[nm] for nm in d.names]
            labels = (torch.as_tensor(f["U"][order], device=dev).to(torch.complex64),
                      torch.as_tensor(f["Z"][order], device=dev).float())

    mcfg = ModelConfig(embed_dim=args.embed_dim, num_layers=args.layers, num_heads=args.heads, n_reps=2,
                       dropout=args.dropout, attention_dropout=args.dropout, predict_residual=True,
                       residual_kappa_scale=args.kscale, residual_z_scale=args.zscale,
                       residual_zero_init=not args.no_zero_init, residual_hyps=args.hyps,
                       residual_hyp_kick=args.hyp_kick)
    model = OTModel(mcfg, mask, args.zhead, input_scale=args.t2_scale).to(dev)
    if args.init_from:
        sd = torch.load(args.init_from, map_location=dev, weights_only=False)
        sd = sd.get("model", sd.get("model_state_dict", sd))
        miss, unexp = model.load_state_dict(sd, strict=False)
        print(f"init from {args.init_from}: missing {len(miss)} unexpected {len(unexp)}")
    unroller = Unroller(args.unroll, args.unroll_eta).to(dev) if args.unroll else None
    params = list(model.parameters()) + (list(unroller.parameters()) if unroller is not None else [])
    nparam = sum(p.numel() for p in model.parameters())
    opt = torch.optim.AdamW(params, lr=args.lr, weight_decay=args.wd)
    steps_per_epoch = math.ceil(len(itr) / args.bs)
    total = args.epochs * steps_per_epoch

    def lr_at(step):
        if step < args.warmup:
            return (step + 1) / args.warmup
        if args.schedule == "constant":
            return 1.0
        p = (step - args.warmup) / max(1, total - args.warmup)
        return 0.5 * (1 + math.cos(math.pi * min(1.0, p)))

    sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_at)
    log({"type": "start", "n_params": nparam, "n_train": len(itr), "n_val": len(iva), "steps": total})
    ev_tr = evaluate(model, d, itr[:141], U_base_all, Z_base_all, mask, args.lam, labels, unroller=unroller)
    ev_va = evaluate(model, d, iva, U_base_all, Z_base_all, mask, args.lam, labels, unroller=unroller)
    log({"type": "eval", "epoch": 0, **{f"train_{k}": v for k, v in ev_tr.items()}, **{f"val_{k}": v for k, v in ev_va.items()}})
    best = float("inf")
    g = torch.Generator(device="cpu").manual_seed(args.seed)
    step = 0
    tgt = None
    pos_in_train = torch.full((d.N,), -1, dtype=torch.long, device=dev)
    pos_in_train[itr] = torch.arange(len(itr), device=dev)
    for ep in range(1, args.epochs + 1):
        t0 = time.time()
        if args.distill_every and ep >= args.distill_start and (ep - args.distill_start) % args.distill_every == 0:
            tA, tS, tD, tres = distill_targets(model, d, itr, U_base_all, Z_base_all, mask, args.lam, args.distill_steps,
                                               args.distill_lr, args.distill_prox, args.kscale)
            tgt = (tA, tS, tD)
            log({"type": "distill", "epoch": ep, "target_resid_median": float(tres.median()), "time_s": time.time() - t0})
        perm = itr[torch.randperm(len(itr), generator=g).to(dev)]
        tot = {"loss": 0.0, "obj": 0.0, "sup": 0.0, "n": 0}
        for s in range(0, len(perm), args.bs):
            b = perm[s:s + args.bs]
            out = model(d.make_batch(b))
            deep = None
            if unroller is not None:
                h = out["residual_hyps"][0]
                pT, traj = unroller([h["dkappa_re"], h["dkappa_im"], h["dz"]], U_base_all[b], Z_base_all[b], mask,
                                    d.t2[b], args.lam, d.znorm_full[b], track=True)
                cands = [assemble(pT, U_base_all[b], Z_base_all[b], mask)]
                deep = torch.stack(traj).mean(0) if traj else None
            else:
                cands = build_params(out, U_base_all[b], Z_base_all[b], mask)
            per = []
            objs, sups = [], []
            for U, Z in cands:
                lk = torch.zeros(len(b), device=dev)
                if args.w_obj:
                    o = per_mol_objective(U, Z, d.t2[b], args.lam, d.znorm_full[b])
                    lk = lk + args.w_obj * o
                    objs.append(o.detach())
                if args.w_sup:
                    sl = sup_loss(U, Z, labels[0][b], labels[1][b], args.w_supz, args.orbit_sup)
                    lk = lk + args.w_sup * sl
                    sups.append(sl.detach())
                per.append(lk)
            P = torch.stack(per)                                    # (K, B)
            lmol = 0.95 * P.min(0).values + 0.05 * P.mean(0) if len(cands) > 1 else P[0]
            if deep is not None and args.unroll_deep:
                lmol = lmol + args.unroll_deep * deep
            if tgt is not None and args.w_distill:
                h = out["residual_hyps"][0]
                j = pos_in_train[b]
                ld = ((h["dkappa_re"] - tgt[0][j]).pow(2).flatten(1).mean(1) + (h["dkappa_im"] - tgt[1][j]).pow(2).flatten(1).mean(1)
                      + (h["dz"] * mask - tgt[2][j]).pow(2).flatten(1).mean(1))
                lmol = lmol + args.w_distill * ld * 100.0
            loss = lmol.mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, args.clip)
            opt.step()
            sched.step()
            step += 1
            tot["loss"] += float(loss) * len(b)
            if objs:
                tot["obj"] += float(torch.stack(objs).min(0).values.sum())
            if sups:
                tot["sup"] += float(torch.stack(sups).min(0).values.sum())
            tot["n"] += len(b)
        rec = {"type": "train", "epoch": ep, "loss": tot["loss"] / tot["n"], "obj": tot["obj"] / tot["n"],
               "sup": tot["sup"] / tot["n"], "lr": sched.get_last_lr()[0], "time": time.time() - t0}
        if ep % args.eval_every == 0 or ep == args.epochs:
            ev_tr = evaluate(model, d, itr[:141], U_base_all, Z_base_all, mask, args.lam, labels, unroller=unroller)
            ev_va = evaluate(model, d, iva, U_base_all, Z_base_all, mask, args.lam, labels, unroller=unroller)
            if unroller is not None:
                rec["etas"] = [round(float(x), 4) for x in torch.exp(unroller.log_eta).flatten()[:9]]
            rec.update({f"train_{k}": v for k, v in ev_tr.items()})
            rec.update({f"val_{k}": v for k, v in ev_va.items()})
            if ev_va["resid_median"] < best:
                best = ev_va["resid_median"]
                torch.save({"model": model.state_dict(), "unroller": unroller.state_dict() if unroller is not None else None,
                            "args": vars(args), "epoch": ep}, out_dir / "best.pt")
        log(rec)
    torch.save({"model": model.state_dict(), "args": vars(args), "epoch": args.epochs}, out_dir / "last.pt")
    log({"type": "done", "best_val_resid_median": best})


if __name__ == "__main__":
    main()
