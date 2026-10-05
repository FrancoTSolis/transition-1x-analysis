#!/usr/bin/env python3
"""GRPO energy fine-tuning of the pretrained chemistry-frame slot model (pretrain.opt_true.slotnet, T = 1).

Policy for molecule i (deterministic network, Gaussian exploration):
    mu_K = SlotNet one-shot generator (real antisymmetric, per rep), mu_Z = zero-init Z head (masked, symmetric)
    a = (a_K, a_Z) ~ N((mu_K, mu_Z), diag(sigma_K^2, sigma_Z^2))          over the free entries only
    U = Phi B expm(a_K)            (chemistry frame B, real sector)
    Z = Z*(U) + a_Z                (Z*: exact VarPro optimum of ffsim's t2 objective for this U)
    r = (E_HF - E(U, Z)) / (E_HF - E_CCSD)       E: exact statevector (norb <= 16) or MPS + SQD
GRPO: group-normalized advantages over G samples per molecule, PPO clipping, optional KL to a periodically
refreshed reference policy (engineering as in pretrain/rl/grpo.py / Transition-State-Generation-Flow).

Usage (from ccsd_amplitudes/):
  python3 -m pretrain.rl.grpo_slot --init runs_ot/slotall_T1/best.pt --train-names <file> --val-names <file> \
      --group-size 8 --batch-mols 4 --steps 50 --n-workers 16 --out rl_runs/slot_dev
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import varpro as V  # noqa: E402
from pretrain.opt_true.eval_slot import load_slot  # noqa: E402
from pretrain.opt_true.train_slot_all import Bucket  # noqa: E402
from pretrain.rl.grpo import _reward_job, _worker_init  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


class ZHead(nn.Module):
    """Zero-init symmetric Z correction read out from the slot model's final pair state."""

    def __init__(self, d, n_reps=2):
        super().__init__()
        self.h = nn.ModuleList([nn.Sequential(nn.Linear(d, d // 2), nn.GELU(), nn.Linear(d // 2, 1)) for _ in range(n_reps)])
        for m in self.h:
            nn.init.zeros_(m[-1].weight)
            nn.init.zeros_(m[-1].bias)

    def forward(self, z):
        zt = z.transpose(1, 2)
        return torch.stack([0.5 * (m(z).squeeze(-1) + m(zt).squeeze(-1)) for m in self.h], 1)   # (B,R,n,n)


class Policy(nn.Module):
    def __init__(self, slot, d):
        super().__init__()
        self.slot = slot
        self.zhead = ZHead(d)

    def forward(self, x, lam):
        """Mean of the final step: (dK real antisym (1,R,n,n), dZ (1,R,n,n)) for a single-molecule batch x.
        Starts from x["Ustart"] (the frame, or the output of k frozen recycles) with recycle state x["prev"]."""
        from pretrain.opt_true.slotnet import features
        pair, node, _, _ = features(x["Ustart"], x["t2"], x["mask"], lam, x["zref"])
        dK, z = self.slot.step(pair, node, x.get("prev"), x.get("t_start", 0), x["chem"])
        return dK.real, self.zhead(z)


@torch.no_grad()
def frozen_prefix(net, x, lam, k):
    """k recycles of the frozen pretrained network: (U_k, recycle state after step k-1)."""
    from pretrain.opt_true.slotnet import features
    U, prev = x["U0"], None
    for t in range(k):
        pair, node, _, _ = features(U, x["t2"], x["mask"], lam, x["zref"])
        dK, prev = net.step(pair, node, prev, t, x["chem"])
        U = U @ torch.linalg.matrix_exp(dK.to(U.dtype))
    return U, prev


def free_idx(n, mask_np):
    iu = np.triu_indices(n, 1)
    pz = np.array([(p, q) for p in range(n) for q in range(p, n) if mask_np[p, q]])
    return iu, (pz[:, 0], pz[:, 1])


def to_flat(dK, dZ, fi):
    (ku, kv), (zu, zv) = fi
    return torch.cat([dK[0][:, ku, kv].reshape(-1), dZ[0][:, zu, zv].reshape(-1)])


def from_flat(a, n, fi, R=2):
    (ku, kv), (zu, zv) = fi
    nk = len(ku)
    K = torch.zeros(R, n, n, dtype=a.dtype, device=a.device)
    K[:, ku, kv] = a[:R * nk].reshape(R, nk)
    K = K - K.transpose(-1, -2)
    Z = torch.zeros(R, n, n, dtype=a.dtype, device=a.device)
    Z[:, zu, zv] = a[R * nk:].reshape(R, len(zu))
    Z = Z + Z.transpose(-1, -2) - torch.diag_embed(torch.diagonal(Z, dim1=-2, dim2=-1))
    return K, Z


class Mol:
    def __init__(self, name, idx, dev, lam, frozen=None, prefix_T=0):
        n, no, nv = idx[name]
        self.name, self.n, self.no = name, n, no
        self.x = Bucket(no, nv, [name]).batch([0], dev)
        d = np.load(ROOT / "rhf_dataset" / f"{name}.npz")
        self.t1 = d["t1"].astype(np.float64)
        self.e_hf, self.e_ccsd = float(d["e_hf"]), float(d["e_ccsd"])
        m = self.x["mask"].cpu().numpy()
        self.fi = free_idx(n, m)
        f = np.load(ROOT / "rhf_frames" / f"{no}_{nv}.npz")
        i = [str(k) for k in f["names"]].index(name)
        sa = f["slot_atom"][i].astype(int)
        bond = f["bond"][i]
        self.bonded_np = (sa[:, None] == sa[None, :]) | bond[sa][:, sa]
        self.nk = len(self.fi[0][0])
        self.lam = lam
        self.x["Ustart"], self.x["prev"], self.x["t_start"] = self.x["U0"], None, 0
        if prefix_T:
            U, prev = frozen_prefix(frozen, self.x, lam, prefix_T)
            self.x["Ustart"], self.x["prev"], self.x["t_start"] = U, prev, prefix_T
        # the frame is stored in float32; ffsim requires unitarity to ~1e-8 -> polar factor in float64
        W, _, Vh = torch.linalg.svd(self.x["Ustart"][0].to(torch.complex128))
        self.U0_64 = W @ Vh

    def realize(self, a):
        """flat action (float64 tensor) -> (U complex128 np (R,n,n), Z float64 np (R,n,n))."""
        K, dZ = from_flat(a, self.n, self.fi)
        U0 = self.U0_64
        U = U0 @ torch.linalg.matrix_exp(K.to(torch.complex128))
        t2 = self.x["t2"].double()
        Zs = V.solve_z(U[None], t2, self.x["mask"].double(), self.lam, self.x["zref"].double())[0]
        Z = (Zs + dZ) * self.x["mask"].double()
        return U.cpu().numpy(), Z.cpu().numpy(), float(D.rel_residual(Z[None], U[None], t2)[0])

    def corr(self, E):
        return (self.e_hf - E) / (self.e_hf - self.e_ccsd)


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--init", required=True, help="slot-model checkpoint (train_slot_all)")
    ap.add_argument("--prefix-T", type=int, default=0, help="frozen recycles before the trainable final step")
    ap.add_argument("--init-policy", default=None,
                    help="resume from a GRPO policy checkpoint (policy_*.pt); --init must be the slot checkpoint it "
                         "was built from (the frozen recycle prefix is taken from --init, i.e. stays pretrained)")
    ap.add_argument("--train-names", required=True)
    ap.add_argument("--val-names", required=True)
    ap.add_argument("--ham-dir", default="rhf_hamiltonians")
    ap.add_argument("--driver-threads", type=int, default=0,
                    help="torch threads for the policy in the driver (0 = leave default); run the job with "
                         "OMP_NUM_THREADS=1 so spawned reward workers start single-threaded")
    ap.add_argument("--tn-chi", type=int, default=64, help="MPS bond dimension for --reward tn")
    ap.add_argument("--tn-basis-cache", default=None, help="shared directory for cached localized bases (--reward tn)")
    ap.add_argument("--tn-stack-mem-gb", type=float, default=1.0)
    ap.add_argument("--tn-cache-items", type=int, default=1, help="TN engines (molecules) cached per reward worker")
    ap.add_argument("--tn-impl", default="current", choices=["current", "v1"],
                    help="v1 = frozen zip-up engine (tn_energy_v1.py) of the Oct-3 n29 run")
    ap.add_argument("--reward", choices=["exact", "mps_sqd", "resid", "queue", "tn"], default="exact",
                    help="resid: 1 - t2 residual (instant; for testing the RL loop mechanics); "
                         "queue: exact energies from GPU workers through pretrain.rl.reward_queue")
    ap.add_argument("--queue-root", default=None, help="shared-FS queue directory for --reward queue")
    ap.add_argument("--queue-timeout", type=float, default=3600.0)
    ap.add_argument("--queue-kind", default="exact", choices=["exact", "tn"],
                    help="kind of the --reward queue tasks: exact (GPU state vector) or tn (served by "
                         "'reward_queue worker --kind tn' workers, whose --chi / --tn-impl set the TN energy)")
    ap.add_argument("--max-bond", type=int, default=32)
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--samples-per-batch", type=int, default=300)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--group-size", type=int, default=8)
    ap.add_argument("--batch-mols", type=int, default=4)
    ap.add_argument("--sigma-k", type=float, default=0.03)
    ap.add_argument("--sigma-z", type=float, default=0.01)
    ap.add_argument("--clip-eps", type=float, default=0.2)
    ap.add_argument("--antithetic", action="store_true", help="mirrored noise pairs (eps, -eps) within each group")
    ap.add_argument("--noise-support", default="full", choices=["full", "bonded"],
                    help="explore only generator entries between same-atom/bonded slots (others keep the mean)")
    ap.add_argument("--ppo-epochs", type=int, default=2)
    ap.add_argument("--kl-coef", type=float, default=0.0)
    ap.add_argument("--ref-update", type=int, default=5)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--zhead-lr", type=float, default=1e-3)
    ap.add_argument("--max-grad-norm", type=float, default=1.0)
    ap.add_argument("--reward-clip", type=float, default=2.0)
    ap.add_argument("--steps", type=int, default=50)
    ap.add_argument("--eval-every", type=int, default=5)
    ap.add_argument("--n-workers", type=int, default=16)
    ap.add_argument("--worker-threads", type=int, default=1)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    torch.manual_seed(args.seed)
    if args.driver_threads:
        torch.set_num_threads(args.driver_threads)
    rng = np.random.default_rng(args.seed)
    dev = args.device
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "args.json").write_text(json.dumps(vars(args), indent=1))
    logf = open(out / "log.jsonl", "a")

    def log(rec):
        rec["time"] = time.time()
        logf.write(json.dumps(rec) + "\n")
        logf.flush()
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    rd = lambda f: [ln.strip() for ln in open(ROOT / f) if ln.strip()]  # noqa: E731
    train_names, val_names = rd(args.train_names), rd(args.val_names)
    slot, a_, _, _ = load_slot(args.init, dev)
    frozen = copy.deepcopy(slot).eval() if args.prefix_T else None
    mols = {n: Mol(n, idx, dev, args.lam, frozen, args.prefix_T) for n in train_names + val_names}
    policy = Policy(slot, a_["d"]).to(dev)
    if args.init_policy:
        ck = torch.load(ROOT / args.init_policy, map_location=dev, weights_only=False)
        assert ck["args"].get("prefix_T", 0) == args.prefix_T, "prefix-T must match the resumed policy"
        policy.load_state_dict(ck["policy"])
        print(f"resumed policy from {args.init_policy} (step {ck.get('step')})", flush=True)
    ref = copy.deepcopy(policy).eval()
    for p in ref.parameters():
        p.requires_grad_(False)
    opt = torch.optim.AdamW([{"params": policy.slot.parameters(), "lr": args.lr},
                             {"params": policy.zhead.parameters(), "lr": args.zhead_lr}], weight_decay=0.0)
    pool = None
    if args.reward in ("exact", "mps_sqd", "tn"):
        pool = get_context("spawn").Pool(args.n_workers, initializer=_worker_init,
                                         initargs=(str(ROOT / args.ham_dir), args.worker_threads))
    if args.reward == "queue":
        from pretrain.rl import reward_queue as RQ
        assert args.queue_root, "--queue-root is required with --reward queue"
    rkw = dict(shots=args.shots, max_bond=args.max_bond, samples_per_batch=args.samples_per_batch)
    if args.reward == "tn":
        rkw = dict(max_bond=args.tn_chi, basis_cache=args.tn_basis_cache, stack_mem=int(args.tn_stack_mem_gb * (1 << 30)),
                   block2_threads=1, tn_cache_items=args.tn_cache_items, impl=args.tn_impl)

    def energies(tasks):
        if args.reward == "resid":            # pseudo-energy with corr fraction = 1 - residual
            out = {}
            for (k, n, U, Z, t1) in tasks:
                m = mols[n]
                r = float(D.rel_residual(torch.as_tensor(Z)[None], torch.as_tensor(U)[None], m.x["t2"].double().cpu())[0])
                out[k] = (m.e_hf - (1.0 - r) * (m.e_hf - m.e_ccsd), {})
            return out
        if args.reward == "queue":
            sub = [{"name": n, "U": U, "Z": Z, "t1": t1, "norb": mols[n].n, "nelec": (mols[n].no, mols[n].no),
                    "kind": args.queue_kind} for (k, n, U, Z, t1) in tasks]
            ids = RQ.submit(args.queue_root, sub)
            # TN energies can run 30 min and their heartbeats come from a child process: be slow to requeue
            res = RQ.collect(args.queue_root, ids, timeout=args.queue_timeout,
                             stale_s=1800.0 if args.queue_kind == "tn" else 300.0)
            return {k: res[i] for (k, *_), i in zip(tasks, ids)}
        jobs = [(k, n, U, Z, t1, "square", args.reward, rkw) for (k, n, U, Z, t1) in tasks]
        return {key: (E, info) for key, E, info in pool.imap_unordered(_reward_job, jobs)}

    def mean_action(pol, m):
        dK, dZ = pol(m.x, args.lam)
        return to_flat(dK.double(), dZ.double(), m.fi)

    def sig_vec(m):
        return torch.cat([torch.full((2 * m.nk,), args.sigma_k), torch.full((2 * len(m.fi[1][0]),), args.sigma_z)]).double().to(dev)

    def explore_mask(m):
        """1 for action entries that are sampled (and enter the likelihood), 0 for deterministic ones."""
        ek = np.ones((2, m.nk), bool)
        if args.noise_support == "bonded":
            (ku, kv), _ = m.fi
            ek[:, ~m.bonded_np[ku, kv]] = False
        return torch.as_tensor(np.concatenate([ek.reshape(-1), np.ones(2 * len(m.fi[1][0]), bool)])).to(dev)

    def sample_noise(n_dim):
        if args.antithetic:
            e = torch.randn(G // 2, n_dim, dtype=torch.float64, device=dev)
            return torch.cat([e, -e])
        return torch.randn(G, n_dim, dtype=torch.float64, device=dev)

    @torch.no_grad()
    def evaluate(names, tag, step):
        policy.eval()
        tasks, res_ = [], {}
        for n in names:
            m = mols[n]
            U, Z, r = m.realize(mean_action(policy, m))
            tasks.append(((n, "pol"), n, U, Z, m.t1))
            res_[n] = r
        res = energies(tasks)
        cf = [mols[n].corr(res[(n, "pol")][0]) for n in names if not math.isnan(res[(n, "pol")][0])]
        log({"type": f"eval_{tag}", "step": step, "n": len(cf), "corr_mean": float(np.mean(cf)),
             "corr_median": float(np.median(cf)), "resid_median": float(np.median(list(res_.values()))),
             "per_mol": {n: round(mols[n].corr(res[(n, 'pol')][0]), 4) for n in names}})
        policy.train()
        return float(np.mean(cf))

    G = args.group_size
    best = evaluate(val_names, "val", 0)
    evaluate(train_names, "train", 0)
    torch.save({"policy": policy.state_dict(), "args": vars(args), "step": 0}, out / "policy_best.pt")
    for step in range(1, args.steps + 1):
        t_step = time.time()
        if args.ref_update > 0 and step % args.ref_update == 0:
            ref.load_state_dict(policy.state_dict())
        batch = list(rng.choice(train_names, size=min(args.batch_mols, len(train_names)), replace=False))
        policy.eval()
        with torch.no_grad():
            mu_old = {n: mean_action(policy, mols[n]) for n in batch}
            acts = {n: mu_old[n][None] + (sig_vec(mols[n]) * explore_mask(mols[n]))[None] * sample_noise(mu_old[n].numel())
                    for n in batch}
            tasks, resid = [], []
            for i, n in enumerate(batch):
                for g in range(G):
                    U, Z, r = mols[n].realize(acts[n][g])
                    tasks.append(((i, g), n, U, Z, mols[n].t1))
                    resid.append(r)
        t_r = time.time()
        res = energies(tasks)
        t_r = time.time() - t_r
        R = torch.full((len(batch), G), float("nan"), dtype=torch.float64)
        for (i, g), (E, info) in res.items():
            if not math.isnan(E):
                R[i, g] = float(np.clip(mols[batch[i]].corr(E), -args.reward_clip, args.reward_clip))
        n_fail = int(torch.isnan(R).sum())
        for i in range(len(batch)):
            ok = ~torch.isnan(R[i])
            R[i][~ok] = R[i][ok].mean() if ok.any() else 0.0
        adv = (R - R.mean(1, keepdim=True)) / R.std(1, keepdim=True, unbiased=False).clamp_min(1e-4)
        adv = adv.to(dev)
        stats = {}
        for ep in range(args.ppo_epochs):
            loss_pg, kl, ratios, act_fr = 0.0, 0.0, [], []
            for i, n in enumerate(batch):
                mu = mean_action(policy, mols[n])
                with torch.no_grad():
                    mu_ref = mean_action(ref, mols[n])
                s2 = sig_vec(mols[n]) ** 2
                em = explore_mask(mols[n]).double()
                a = acts[n]
                lr_ = ((((a - mu_old[n][None]) ** 2 - (a - mu[None]) ** 2) / (2 * s2[None])) * em[None]).sum(-1)
                ratio = torch.exp(lr_.clamp(-20, 20))
                A = adv[i]
                surr = torch.minimum(ratio * A, ratio.clamp(1 - args.clip_eps, 1 + args.clip_eps) * A)
                act_fr.append((((A >= 0) & (ratio <= 1 + args.clip_eps)) | ((A < 0) & (ratio >= 1 - args.clip_eps))).float())
                loss_pg = loss_pg - surr.mean()
                kl = kl + ((((mu - mu_ref) ** 2) / (2 * s2)) * em).sum()
                ratios.append(ratio.detach())
            loss = (loss_pg + args.kl_coef * kl) / len(batch)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            gn = torch.nn.utils.clip_grad_norm_(policy.parameters(), args.max_grad_norm)
            opt.step()
            rt = torch.cat(ratios)
            stats = dict(loss=float(loss), kl=float(kl / len(batch)), grad_norm=float(gn), ratio_mean=float(rt.mean()),
                         ratio_max=float(rt.max()), clip_active=float(torch.cat(act_fr).mean()))
        log(dict(type="train", step=step, reward_mean=float(R.mean()), reward_max=float(R.max()),
                 reward_std_in_group=float(R.std(1).mean()), resid_mean=float(np.mean(resid)), n_fail=n_fail,
                 t_reward=t_r, t_step=time.time() - t_step, **stats))
        if step % args.eval_every == 0 or step == args.steps:
            v = evaluate(val_names, "val", step)
            torch.save({"policy": policy.state_dict(), "args": vars(args), "step": step}, out / "policy_last.pt")
            if v > best:
                best = v
                torch.save({"policy": policy.state_dict(), "args": vars(args), "step": step}, out / "policy_best.pt")
    evaluate(train_names, "train", args.steps)
    if pool is not None:
        pool.close()
        pool.join()
    log({"type": "done", "best_val_corr_mean": best})


if __name__ == "__main__":
    main()
