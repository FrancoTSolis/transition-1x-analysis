#!/usr/bin/env python3
"""GRPO fine-tuning of the residual LUCJ policy with an energy reward.

Engineering borrowed from Transition-State-Generation-Flow/InterpolationFlow/
pl_modules/model.py::compute_rl_loss (Gaussian policy with fixed sigma,
group-normalized advantages, PPO clipping, optional KL to a periodically
refreshed reference policy, optional anchor to the pretraining objective, and
the same diagnostic logging: reward, ratio, clip-active fraction, adv sign
fractions, KL).

Rollout for one molecule i (input t2_i):
    mu_i      = policy(t2_i)                      flat residual (dkappa, dZ)
    a_ig      = mu_old_i + sigma * eps_ig,   g = 1..G      (group)
    (U, Z)_ig = canonical_exact_init(t2_i) + residual(a_ig)
    E_ig      = energy(U, Z)      [exact statevector | MPS + SQD], CPU pool
    r_ig      = (E_HF - E_ig) / (E_HF - E_CCSD)   correlation energy recovered
    A_ig      = (r_ig - mean_g r) / (std_g r + eps)
    ratio_ig  = exp( (||a - mu_old||^2 - ||a - mu||^2) / (2 sigma^2) )
    L = -mean min(ratio A, clip(ratio) A) + kl_coef * ||mu - mu_ref||^2/(2 sigma^2)
        + anchor_coef * mean ||mu||^2

Usage (train venv, from ccsd_amplitudes/):
    python -m pretrain.rl.grpo --names-file gauge_study/names_small_norb_le18.txt \
        --reward exact --group-size 8 --batch-mols 4 --sigma 0.05 --steps 100 \
        --init-backbone checkpoints_invariant/best.pt --out-dir rl_runs/dev
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")   # inherited by spawned reward workers
import numpy as np  # noqa: E402
import torch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gauge_study.compressed_canonical import canonical_exact_init, z_mask  # noqa: E402
from pretrain.dataset import CCSDAmplitudeDataset  # noqa: E402
from pretrain.model import ModelConfig  # noqa: E402
from pretrain.rl.policy import (PolicyConfig, ResidualPolicy, apply_residual,  # noqa: E402
                                free_index, gather_flat, unflatten_np)


# ------------------------------------------------------------ reward pool

_HAM_CACHE: dict = {}


def _worker_init(ham_dir: str, n_threads: int):
    # RAYON/NUMBA too: ffsim's compiled kernels use their own thread pool (measured on scai1: unpinned workers
    # take ~22 cores each, 31 s/energy vs 199 s single-threaded -> 3.5x worse core-efficiency, and oversubscription)
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[_v] = str(n_threads)
    try:
        import torch as _t
        _t.set_num_threads(n_threads)
    except Exception:  # noqa: BLE001
        pass
    global _HAM_DIR
    _HAM_DIR = ham_dir


def _get_ham(name: str):
    if name not in _HAM_CACHE:
        from pretrain.rl.hamiltonian import load_hamiltonian
        _HAM_CACHE[name] = load_hamiltonian(_HAM_DIR, name)
    return _HAM_CACHE[name]


_TN_CACHE: dict = {}


def _tn_energy(name, U, Z, t1, kw):
    """MPS (tensor-network) LUCJ energy, pretrain.rl.tn_energy; engines cached per (molecule, chi) in this process.
    kw["impl"] == "v1" selects the frozen zip-up engine (tn_energy_v1.py) that the Oct-3 n29 RL run used."""
    import tempfile
    if kw.get("impl", "current") == "v1":
        from pretrain.rl.tn_energy_v1 import LUCJEnergyTN
    else:
        from pretrain.rl.tn_energy import LUCJEnergyTN
    chi = int(kw.get("max_bond", 64))
    key = (name, chi)
    if key not in _TN_CACHE:
        while len(_TN_CACHE) >= int(kw.get("tn_cache_items", 1)):
            _TN_CACHE.pop(next(iter(_TN_CACHE)))
            import gc
            gc.collect()
        d = np.load(Path(_HAM_DIR) / f"{name}.npz")
        scratch = tempfile.mkdtemp(prefix=f"tn_{os.getpid()}_", dir=os.environ.get("TMPDIR", None))
        _TN_CACHE[key] = LUCJEnergyTN(d["one_body"], d["two_body"], float(d["constant"]), int(d["norb"]),
                                     (int(d["nelec_a"]), int(d["nelec_b"])), max_bond=chi, device="cpu", name=name,
                                     block2_threads=int(kw.get("block2_threads", 1)), scratch=scratch,
                                     stack_mem=int(kw.get("stack_mem", 1 << 30)), basis_cache=kw.get("basis_cache"))
    E, info = _TN_CACHE[key].energy(U, Z, t1=t1)
    return float(E), {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                      for k, v in info.items() if k in ("discarded_sum", "t_total", "max_bond")}


def _reward_job(task):
    """task: (key, name, U, Z, t1, connectivity, kind, kwargs) -> (key, E, info)."""
    key, name, U, Z, t1, conn, kind, kw = task
    if kind == "tn":
        try:
            t0 = time.time()
            E, info = _tn_energy(name, U, Z, t1, kw)
            info["t"] = time.time() - t0
            return key, E, info
        except Exception as e:  # noqa: BLE001
            return key, float("nan"), {"error": f"{type(e).__name__}: {e}"}
    try:
        from pretrain.rl.energy import exact_energy, make_ucj_op, mps_sqd_energy
        ham, norb, nelec, e_hf, e_ccsd = _get_ham(name)
        op = make_ucj_op(Z, U, conn, t1=t1)
        t0 = time.time()
        if kind == "exact":
            E = exact_energy(ham, norb, nelec, op)
            info = {}
        else:
            E, info = mps_sqd_energy(ham, norb, nelec, op, **kw)
        info["t"] = time.time() - t0
        return key, E, info
    except Exception as e:  # noqa: BLE001
        return key, float("nan"), {"error": f"{type(e).__name__}: {e}"}


# ------------------------------------------------------------------ utils

def load_names(path: str, ham_dir: str, data_dir: str) -> list[str]:
    names = [ln.strip() for ln in open(path) if ln.strip()]
    return [n for n in names if (Path(ham_dir) / f"{n}.npz").exists()
            and (Path(data_dir) / f"{n}.npz").exists()]


def mol_context(ds: CCSDAmplitudeDataset, name: str, connectivity: str):
    d = np.load(ds.data_dir / f"{name}.npz")
    t2 = d["t2"].astype(np.float64)
    init = canonical_exact_init(t2)
    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    m = z_mask(connectivity, norb)
    return dict(name=name, t2=t2, t1=d["t1"].astype(np.float64), nocc=nocc, nvirt=nvirt,
                norb=norb, U_init=init.U, Z_init=init.Z * m[None], mask=m,
                fidx=free_index(norb, m), e_hf=float(d["e_hf"]), e_ccsd=float(d["e_ccsd"]))


def flat_mu(out, batch, ctxs, device):
    """Per-sample flat mean vectors (list of 1-D tensors, with grad)."""
    max_nocc = batch["max_nocc"]
    mus = []
    for i, c in enumerate(ctxs):
        idx = torch.cat([torch.arange(c["nocc"], device=device),
                         max_nocc + torch.arange(c["nvirt"], device=device)])
        mus.append(gather_flat(out, i, idx, c["fidx"]))
    return mus


# ------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--ham-dir", default="rhf_hamiltonians")
    ap.add_argument("--names-file", required=True)
    ap.add_argument("--val-frac", type=float, default=0.2)
    ap.add_argument("--connectivity", default="square")
    ap.add_argument("--reward", choices=["exact", "mps_sqd"], default="exact")
    ap.add_argument("--max-bond", type=int, default=64)
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--samples-per-batch", type=int, default=300)
    ap.add_argument("--reward-clip", type=float, default=2.0,
                    help="clip |reward| (corr. energy fraction) for stability")
    # policy / init
    ap.add_argument("--init-backbone", default=None, help="pretrained PretrainingModel ckpt")
    ap.add_argument("--init-policy", default=None, help="resume a full policy ckpt")
    ap.add_argument("--embed-dim", type=int, default=192)
    ap.add_argument("--num-layers", type=int, default=6)
    ap.add_argument("--num-heads", type=int, default=8)
    ap.add_argument("--kappa-scale", type=float, default=1.0)
    ap.add_argument("--z-scale", type=float, default=1.0)
    ap.add_argument("--freeze-backbone", action="store_true")
    # GRPO
    ap.add_argument("--group-size", type=int, default=8)
    ap.add_argument("--batch-mols", type=int, default=4)
    ap.add_argument("--sigma", type=float, default=0.05)
    ap.add_argument("--sample-from", choices=["current", "ref"], default="current")
    ap.add_argument("--ppo-epochs", type=int, default=1)
    ap.add_argument("--clip-eps", type=float, default=0.2)
    ap.add_argument("--kl-coef", type=float, default=0.0)
    ap.add_argument("--ref-update", type=int, default=5, help="refresh ref policy every N steps (0=never)")
    ap.add_argument("--anchor-coef", type=float, default=0.0, help="coef of mean ||mu||^2")
    ap.add_argument("--lr", type=float, default=1e-5)
    ap.add_argument("--max-grad-norm", type=float, default=1.0)
    ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--eval-every", type=int, default=5)
    ap.add_argument("--n-workers", type=int, default=32)
    ap.add_argument("--worker-threads", type=int, default=1)
    ap.add_argument("--out-dir", default="rl_runs/dev")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    log_f = open(out_dir / "log.jsonl", "a")
    (out_dir / "args.json").write_text(json.dumps(vars(args), indent=1))

    def log(rec):
        rec["time"] = time.time()
        log_f.write(json.dumps(rec) + "\n"); log_f.flush()
        print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

    # ---- data
    names = load_names(args.names_file, args.ham_dir, args.data_dir)
    rng.shuffle(names)
    n_val = max(1, int(len(names) * args.val_frac))
    val_names, train_names = names[:n_val], names[n_val:]
    ds = CCSDAmplitudeDataset(args.data_dir)
    name_to_idx = {n: i for i, n in enumerate(ds.names)}
    ctx_cache: dict[str, dict] = {}

    def ctx(n):
        if n not in ctx_cache:
            ctx_cache[n] = mol_context(ds, n, args.connectivity)
        return ctx_cache[n]

    print(f"train {len(train_names)} / val {len(val_names)} molecules; reward={args.reward}; device={device}")

    # ---- policy
    # no dropout anywhere: the rollout mean mu_old (eval mode) and the update
    # mean mu must come from the same deterministic network, otherwise the
    # Gaussian log-ratio at sigma~0.02 explodes (observed: ratio -> 0, KL 246).
    mcfg = ModelConfig(embed_dim=args.embed_dim, num_layers=args.num_layers,
                       num_heads=args.num_heads, n_reps=2, dropout=0.0,
                       attention_dropout=0.0, predict_invariant=True)
    pcfg = PolicyConfig(kappa_scale=args.kappa_scale, z_scale=args.z_scale,
                        zero_init=True, connectivity=args.connectivity)
    policy = ResidualPolicy(mcfg, pcfg).to(device)
    if args.init_policy:
        policy.load_state_dict(torch.load(args.init_policy, map_location=device)["policy"])
    elif args.init_backbone:
        sd = torch.load(args.init_backbone, map_location=device, weights_only=False)["model_state_dict"]
        miss, unexp = policy.load_backbone(sd)
        print(f"loaded backbone from {args.init_backbone} (missing {len(miss)}, unexpected {len(unexp)})")
        # a compressed_recon / compressed_sup checkpoint carries residual heads:
        # start the policy mean at the amortized-optimizer solution.
        res_sd = {k[len("decode_heads.residual."):]: v for k, v in sd.items()
                  if k.startswith("decode_heads.residual.")}
        if res_sd:
            policy.heads.load_state_dict(res_sd)
            print(f"loaded residual heads from checkpoint ({len(res_sd)} tensors)")
    import copy
    ref = copy.deepcopy(policy).eval()
    for p in ref.parameters():
        p.requires_grad_(False)
    params = list(policy.heads.parameters()) if args.freeze_backbone else list(policy.parameters())
    opt = torch.optim.AdamW(params, lr=args.lr, weight_decay=0.0)

    # ---- reward pool
    ctxmp = get_context("spawn")
    pool = ctxmp.Pool(args.n_workers, initializer=_worker_init,
                      initargs=(args.ham_dir, args.worker_threads))
    rkw = dict(shots=args.shots, max_bond=args.max_bond, samples_per_batch=args.samples_per_batch)

    def energies(tasks):
        """tasks: list of (key, name, U, Z, t1). Returns dict key -> (E, info)."""
        jobs = [(k, n, U, Z, t1, args.connectivity, args.reward, rkw) for (k, n, U, Z, t1) in tasks]
        res = {}
        for key, E, info in pool.imap_unordered(_reward_job, jobs):
            res[key] = (E, info)
        return res

    def make_batch(mol_names):
        samples = [ds[name_to_idx[n]] for n in mol_names]
        b = CCSDAmplitudeDataset.collate_fn(samples)
        return {k: v.to(device) if isinstance(v, torch.Tensor) else v for k, v in b.items()}

    def corr_frac(E, c):
        return (c["e_hf"] - E) / (c["e_hf"] - c["e_ccsd"])

    @torch.no_grad()
    def evaluate(mol_names, tag, step):
        policy.eval()
        tasks = []
        for s in range(0, len(mol_names), 16):
            chunk = mol_names[s:s + 16]
            b = make_batch(chunk)
            out = policy(b)
            mus = flat_mu(out, b, [ctx(n) for n in chunk], device)
            for n, mu in zip(chunk, mus):
                c = ctx(n)
                dk, dZ = unflatten_np(mu.cpu().numpy().astype(np.float64), c["norb"], c["fidx"])
                U, Z = apply_residual(c["U_init"], c["Z_init"], dk, dZ)
                tasks.append(((n, "pol"), n, U, Z, c["t1"]))
                if step == 0:
                    tasks.append(((n, "init"), n, c["U_init"], c["Z_init"], c["t1"]))
        res = energies(tasks)
        fr_pol = [corr_frac(res[(n, "pol")][0], ctx(n)) for n in mol_names if not math.isnan(res[(n, "pol")][0])]
        rec = {"type": f"eval_{tag}", "step": step, "n": len(fr_pol),
               "corr_frac_mean": float(np.mean(fr_pol)), "corr_frac_median": float(np.median(fr_pol))}
        if step == 0:
            fr_init = [corr_frac(res[(n, "init")][0], ctx(n)) for n in mol_names]
            rec["init_corr_frac_mean"] = float(np.mean(fr_init))
            rec["init_corr_frac_median"] = float(np.median(fr_init))
        log(rec)
        policy.train()

    # ---- training loop
    evaluate(val_names, "val", 0)
    evaluate(train_names[:min(len(train_names), 32)], "train", 0)
    G, sig = args.group_size, args.sigma
    for step in range(1, args.steps + 1):
        t_step = time.time()
        if args.ref_update > 0 and step % args.ref_update == 0:
            ref.load_state_dict(policy.state_dict())
        mol_names = list(rng.choice(train_names, size=min(args.batch_mols, len(train_names)), replace=False))
        ctxs = [ctx(n) for n in mol_names]
        b = make_batch(mol_names)

        # ---- rollout (sample actions from old policy)
        policy.eval()
        with torch.no_grad():
            out_old = (policy if args.sample_from == "current" else ref)(b)
            mu_old = [m.detach() for m in flat_mu(out_old, b, ctxs, device)]
            actions = [mo[None, :] + sig * torch.randn(G, mo.numel(), device=device) for mo in mu_old]
        policy.train()

        # ---- rewards
        t_r = time.time()
        tasks = []
        for i, c in enumerate(ctxs):
            A = actions[i].cpu().numpy().astype(np.float64)
            for g in range(G):
                dk, dZ = unflatten_np(A[g], c["norb"], c["fidx"])
                U, Z = apply_residual(c["U_init"], c["Z_init"], dk, dZ)
                tasks.append(((i, g), c["name"], U, Z, c["t1"]))
        res = energies(tasks)
        t_r = time.time() - t_r
        R = torch.zeros(len(ctxs), G, device=device)
        n_fail = 0
        for (i, g), (E, info) in res.items():
            if math.isnan(E):
                n_fail += 1
                R[i, g] = float("nan")
            else:
                R[i, g] = float(np.clip(corr_frac(E, ctxs[i]), -args.reward_clip, args.reward_clip))
        # failed evaluations get the group mean (zero advantage)
        for i in range(len(ctxs)):
            row = R[i]; ok = ~torch.isnan(row)
            if ok.sum() == 0:
                R[i] = 0.0
            else:
                R[i][~ok] = row[ok].mean()
        adv = R - R.mean(1, keepdim=True)
        adv = adv / adv.std(1, keepdim=True, unbiased=False).clamp_min(1e-4)

        # ---- policy update(s)  (eval mode = deterministic; grads still flow)
        policy.eval()
        stats = {}
        for ep in range(args.ppo_epochs):
            out = policy(b)
            mus = flat_mu(out, b, ctxs, device)
            with torch.no_grad():
                mus_ref = flat_mu(ref(b), b, ctxs, device)
            loss_pg = 0.0; kl = 0.0; anchor = 0.0
            ratios, active = [], []
            for i in range(len(ctxs)):
                a = actions[i]                                              # (G, P)
                lr_ = ((a - mu_old[i][None]).pow(2).sum(-1) - (a - mus[i][None]).pow(2).sum(-1)) / (2 * sig ** 2)
                ratio = torch.exp(lr_.clamp(-20, 20))                        # (G,)
                A = adv[i].detach()
                surr = ratio * A
                if args.clip_eps is not None and args.clip_eps > 0:
                    surr_c = ratio.clamp(1 - args.clip_eps, 1 + args.clip_eps) * A
                    act = ((A >= 0) & (ratio <= 1 + args.clip_eps)) | ((A < 0) & (ratio >= 1 - args.clip_eps))
                    surr = torch.minimum(surr, surr_c)
                    active.append(act.float())
                loss_pg = loss_pg - surr.mean()
                kl = kl + (mus[i] - mus_ref[i]).pow(2).sum() / (2 * sig ** 2)
                anchor = anchor + mus[i].pow(2).mean()
                ratios.append(ratio.detach())
            nB = len(ctxs)
            loss = loss_pg / nB + args.kl_coef * kl / nB + args.anchor_coef * anchor / nB
            opt.zero_grad(set_to_none=True)
            loss.backward()
            gn = torch.nn.utils.clip_grad_norm_(params, args.max_grad_norm)
            opt.step()
            ratios_t = torch.cat(ratios)
            stats = dict(loss=float(loss), loss_pg=float(loss_pg / nB), kl=float(kl / nB),
                         anchor=float(anchor / nB), grad_norm=float(gn),
                         ratio_mean=float(ratios_t.mean()), ratio_max=float(ratios_t.max()),
                         clip_active_frac=float(torch.cat(active).mean()) if active else 1.0)
        mu_norm = float(np.mean([m.detach().norm().item() / math.sqrt(m.numel()) for m in mus]))
        rec = dict(type="train", step=step, reward_mean=float(R.mean()), reward_std_in_group=float(R.std(1).mean()),
                   reward_max=float(R.max()), adv_pos_frac=float((adv > 0).float().mean()),
                   n_fail=n_fail, mu_rms=mu_norm, t_reward=t_r, t_step=time.time() - t_step, **stats)
        log(rec)
        if step % args.eval_every == 0:
            evaluate(val_names, "val", step)
            torch.save({"policy": policy.state_dict(), "args": vars(args), "step": step},
                       out_dir / "policy_last.pt")
    pool.close(); pool.join()


if __name__ == "__main__":
    main()
