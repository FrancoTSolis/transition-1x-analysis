#!/usr/bin/env python3
"""Exact LUCJ energies of pretrained optimize=True models on the small molecules (norb <= 16).

For each molecule (any shape; the models are size-agnostic pair-token networks) this builds
  init      canonical exact-DF init (U0, masked Z0)                         [what ffsim starts from]
  net       network output (canonical frame, right-multiplied generator, pair Z)
  net+unroll  network + its learned unrolled steps (if the checkpoint has an unroller)
  net+adamS   network output refined by S Adam steps on ffsim's complex objective
  init+adamS  canonical init refined by S Adam steps (a GPU stand-in for ffsim's optimize=True)
and evaluates the exact statevector energy of the n_reps=2 square-connectivity LUCJ state (with the t1
final orbital rotation), reported as % of the CCSD correlation energy recovered.

Usage (train venv; CPU pool for energies):
  python3 -m pretrain.opt_true.eval_energy --ckpt runs_ot/full_unroll8/best.pt --tag unroll8 \
      --names-file gauge_study/names_small_norb_le16.txt --adam-steps 500 --n-procs 24
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gauge_study.compressed_canonical import canonical_exact_init  # noqa: E402
from pretrain.model import ModelConfig  # noqa: E402
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import train as T  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def mol_batch(name, dev):
    d = np.load(ROOT / "rhf_dataset" / f"{name}.npz")
    t2 = d["t2"].astype(np.float64)
    nocc, _, nvirt, _ = t2.shape
    n = nocc + nvirt
    init = canonical_exact_init(t2)
    m = D.square_mask(n, dev)
    b = {"t2": torch.as_tensor(t2, device=dev).float()[None],
         "noccs": torch.tensor([nocc], device=dev), "nvirts": torch.tensor([nvirt], device=dev),
         "norbs": torch.tensor([n], device=dev), "nocc": [nocc], "nvirt": [nvirt], "n_reps": [2],
         "max_nocc": nocc, "max_nvirt": nvirt}
    U0 = torch.as_tensor(init.U, device=dev).to(torch.complex64)[None]
    Z0 = (torch.as_tensor(init.Z, device=dev).float() * m)[None]
    zref = torch.tensor([init.znorm_full], device=dev).float()
    return b, U0, Z0, m, zref, d["t1"].astype(np.float64), float(d["e_hf"]), float(d["e_ccsd"])


def load_model(ckpt, dev):
    sd = torch.load(ckpt, map_location=dev, weights_only=False)
    a = sd["args"]
    mcfg = ModelConfig(embed_dim=a["embed_dim"], num_layers=a["layers"], num_heads=a["heads"], n_reps=2, dropout=0.0,
                       attention_dropout=0.0, predict_residual=True, residual_kappa_scale=a["kscale"],
                       residual_z_scale=a["zscale"], residual_zero_init=True, residual_hyps=a["hyps"],
                       residual_hyp_kick=a["hyp_kick"])
    mask = D.square_mask(29, dev)
    model = T.OTModel(mcfg, mask, a.get("zhead", "pair"), input_scale=a.get("t2_scale", 1.0)).to(dev)
    model.load_state_dict(sd["model"])
    model.eval()
    unroller = None
    if sd.get("unroller"):
        unroller = T.Unroller(a["unroll"], a["unroll_eta"]).to(dev)
        unroller.load_state_dict(sd["unroller"])
        unroller.eval()
    T.PARAM_MODE = a.get("param", "right")
    assert a.get("frame", "canonical") == "canonical" and a.get("zhead", "pair") == "pair", "only canonical/pair models"
    return model, unroller, a


def energy_job(task):
    key, name, U, Z, t1 = task
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from pretrain.rl.energy import exact_energy, make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
    # float32 predictions are unitary only to ~1e-6; ffsim checks tightly -> nearest unitary (polar, float64)
    W, _, Vh = np.linalg.svd(U.astype(np.complex128))
    U = W @ Vh
    t0 = time.time()
    E = exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1=t1))
    return key, E, (e_hf - E) / (e_hf - e_ccsd), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", action="append", default=[], help="opt_true checkpoint(s); repeatable")
    ap.add_argument("--tag", action="append", default=[], help="name per --ckpt")
    ap.add_argument("--names-file", default="gauge_study/names_small_norb_le16.txt")
    ap.add_argument("--adam-steps", type=int, default=500)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--n-procs", type=int, default=24)
    ap.add_argument("--out", default="pretrain/opt_true/results/energy_small.json")
    ap.add_argument("--skip-init", action="store_true")
    args = ap.parse_args()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    names = [ln.strip() for ln in open(ROOT / args.names_file) if ln.strip()]
    names = [n for n in names if (ROOT / "rhf_hamiltonians" / f"{n}.npz").exists()]
    tasks, meta = [], {}
    t0 = time.time()
    for name in names:
        b, U0, Z0, m, zref, t1, e_hf, e_ccsd = mol_batch(name, dev)
        t2 = b["t2"]
        resid = {}
        cands = {}
        if not args.skip_init:
            cands["init"] = (U0, Z0)
            U, Z, _ = D.refine(t2, U0, Z0, m, lam=args.lam, znorm_ref=zref, steps=args.adam_steps, lr=0.01)
            cands[f"init+adam{args.adam_steps}"] = (U, Z)
        for ck, tag in zip(args.ckpt, args.tag):
            model, unroller, a = load_model(ck, dev)
            with torch.no_grad():
                h = model(b)["residual_hyps"][0]
            p0 = [h["dkappa_re"], h["dkappa_im"], h["dz"]]
            U, Z = T.assemble(p0, U0, Z0, m)
            cands[f"{tag}"] = (U.detach(), Z.detach())
            if unroller is not None:
                with torch.enable_grad():
                    pp = [x.detach().requires_grad_(True) for x in p0]
                    pT, _ = unroller(pp, U0, Z0, m, t2, args.lam, zref)
                U2, Z2 = T.assemble([x.detach() for x in pT], U0, Z0, m)
                cands[f"{tag}+unroll"] = (U2.detach(), Z2.detach())
                base = [x.detach() for x in pT]
            else:
                base = [x.detach() for x in p0]
            with torch.enable_grad():
                U3, Z3, _ = D.refine(t2, U0, Z0, m, lam=args.lam, znorm_ref=zref, steps=args.adam_steps, lr=0.01,
                                     A0=base[0], S0=base[1], D0=base[2])
            cands[f"{tag}+adam{args.adam_steps}"] = (U3, Z3)
        for k, (U, Z) in cands.items():
            resid[k] = float(D.rel_residual(Z.double(), U.to(torch.complex128), t2.double())[0])
            tasks.append(((name, k), name, U[0].detach().cpu().numpy().astype(np.complex128),
                          Z[0].detach().cpu().numpy().astype(np.float64), t1))
        meta[name] = {"resid": resid, "e_hf": e_hf, "e_ccsd": e_ccsd, "norb": int(b["norbs"][0])}
    print(f"{len(tasks)} energies for {len(names)} molecules ({time.time()-t0:.0f}s to build)", flush=True)
    ctx = get_context("spawn")
    res = {}
    with ctx.Pool(args.n_procs) as pool:
        for (name, k), E, cf, dt in pool.imap_unordered(energy_job, tasks):
            res.setdefault(name, {})[k] = {"E": E, "corr_frac": cf, "t": dt}
            print(f"  {name:22s} {k:28s} corr% {100*cf:7.1f}  resid {meta[name]['resid'][k]:.3f}  ({dt:.0f}s)", flush=True)
    keys = list(next(iter(res.values())).keys())
    print("\n=== median over molecules: % CCSD correlation energy recovered (resid) ===")
    summ = {}
    for k in sorted(keys):
        cf = np.array([100 * res[n][k]["corr_frac"] for n in res if k in res[n]])
        rs = np.array([meta[n]["resid"][k] for n in res if k in res[n]])
        summ[k] = {"median_corr_pct": float(np.median(cf)), "mean_corr_pct": float(np.mean(cf)),
                   "min_corr_pct": float(cf.min()), "median_resid": float(np.median(rs)), "n": int(len(cf))}
        print(f"  {k:30s} median {np.median(cf):6.1f}  mean {np.mean(cf):6.1f}  min {cf.min():7.1f}  resid {np.median(rs):.3f}  (n={len(cf)})")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump({"summary": summ, "per_molecule": res, "meta": meta, "args": vars(args)}, open(args.out, "w"), indent=1)
    print(f"-> {args.out}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
