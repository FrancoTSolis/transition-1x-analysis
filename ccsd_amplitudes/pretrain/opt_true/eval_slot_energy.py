#!/usr/bin/env python3
"""Exact LUCJ energies (statevector, ffsim) on the 30 held-out small molecules (norb 15-16) for:
  init            canonical exact-DF init, masked Z (what ffsim's optimize=True starts from)
  label           ffsim optimize=True from the canonical init (rhf_targets_compressed_small, L-BFGS 500)
  frame           chemistry frame (A = 0) + exact VarPro Z
  frame+adamN     chemistry frame + N in-frame Adam steps on real generators (full support)
  <tag>_t<t>      slot model after t recycles (one network call per recycle)
  <tag>_t<t>+adamM  that output + M in-frame Adam steps
Reports % of the CCSD correlation energy recovered, (E_HF - E) / (E_HF - E_CCSD), and the t2 residual.

Usage: python3 -m pretrain.opt_true.eval_slot_energy --ckpt runs_ot/slotall_T1/best.pt --tag all1 [--ckpt ... --tag ...]
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
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true import frame_eval as FE  # noqa: E402
from pretrain.opt_true import varpro as V  # noqa: E402
from pretrain.opt_true.eval_slot import load_slot  # noqa: E402
from pretrain.opt_true.train_slot_all import Bucket  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def energy_job(task):
    key, name, U, Z, t1 = task
    nt = os.environ.get("ENERGY_THREADS", "1")   # per-worker threads (large norb: few workers, many threads)
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[v] = nt                    # ffsim's own thread pool ignores OMP_NUM_THREADS
    from pretrain.rl.energy import exact_energy, make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
    W, _, Vh = np.linalg.svd(U.astype(np.complex128))
    U = W @ Vh
    t0 = time.time()
    E = exact_energy(ham, norb, nelec, make_ucj_op(Z, U, "square", t1=t1))
    return key, E, (e_hf - E) / (e_hf - e_ccsd), time.time() - t0


def in_frame_adam(U, t2, mask, lam, zref, steps, lr=0.02, nocc=None):
    """Continue from a real-sector U = Phi O by Adam on real generators (full support), exact VarPro Z."""
    n = U.shape[-1]
    ph = FE.real_sector_phases(n, nocc, U.device)
    O = (ph.conj()[None, :, :, None] * U.to(torch.complex128)).real          # Phi^-1 U, real orthogonal
    A, _ = FE.optimize(O, ph, t2, mask, lam, zref, [torch.ones_like(O)], [steps], lr=lr)
    Uo = FE.u_from_real(O, A, ph)
    return Uo, V.solve_z(Uo, t2, mask, lam, zref)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", action="append", default=[])
    ap.add_argument("--tag", action="append", default=[])
    ap.add_argument("--T", type=int, nargs="*", default=None, help="recycle counts to evaluate (default: 1 and train T)")
    ap.add_argument("--adam", type=int, default=300, help="in-frame Adam steps for frame+adam")
    ap.add_argument("--net-adam", type=int, default=100, help="in-frame Adam steps after the network")
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--n-procs", type=int, default=12)
    ap.add_argument("--names-file", default="gauge_study/names_small_norb_le16.txt")
    ap.add_argument("--labels-dir", default="rhf_targets_compressed_small")
    ap.add_argument("--out", default="pretrain/opt_true/results/energy_slot_small.json")
    ap.add_argument("--skip-energy", action="store_true")
    ap.add_argument("--no-baselines", action="store_true", help="only the network candidates")
    ap.add_argument("--dump", default=None, help="build candidates (GPU) and pickle the energy tasks; no energies")
    ap.add_argument("--from-dump", action="append", default=[], help="compute energies for pickled tasks (CPU only)")
    ap.add_argument("--skip-keys", nargs="*", default=[], help="candidate keys not to compute (e.g. label)")
    args = ap.parse_args()
    import pickle
    if args.from_dump:
        tasks, meta = [], {}
        for f in args.from_dump:
            dd = pickle.load(open(f, "rb"))
            tasks += dd["tasks"]
            for n, m in dd["meta"].items():
                meta.setdefault(n, {"resid": {}, "norb": m["norb"]})["resid"].update(m["resid"])
        names = list(meta)
        keys = list(dict.fromkeys(k for n in meta for k in meta[n]["resid"]))
        tasks = list({t[0]: t for t in tasks}.values())                   # dedupe (label appears in several dumps)
        tasks = [t for t in tasks if t[0][1] not in args.skip_keys]
        keys = [k for k in keys if k not in args.skip_keys]
        return run_energies(args, tasks, meta, names, keys, time.time())
    dev = "cuda"
    names = [ln.strip() for ln in open(ROOT / args.names_file) if ln.strip()]
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    models = [(tag, *load_slot(ck, dev)) for ck, tag in zip(args.ckpt, args.tag)]
    tasks, meta = [], {}
    t0 = time.time()
    for name in names:
        n, no, nv = idx[name]
        bk = Bucket(no, nv, [name])
        x = bk.batch([0], dev)
        t2 = x["t2"].double()
        zr = x["zref"].double()
        mask = D.square_mask(n, dev, torch.float64)
        dd = np.load(ROOT / "rhf_dataset" / f"{name}.npz")
        t1 = dd["t1"].astype(np.float64)
        cands = {}
        U0 = x["U0"].to(torch.complex128)
        ini = np.load(ROOT / args.labels_dir / "init" / f"{name}.npz")
        cands["init"] = (torch.as_tensor(ini["U_re"] + 1j * ini["U_im"], device=dev)[None],
                         torch.as_tensor(ini["Z"], device=dev)[None].double() * mask)
        lab = np.load(ROOT / args.labels_dir / "square_reg0.005" / f"{name}.npz")
        cands["label"] = (torch.as_tensor(lab["U_re"] + 1j * lab["U_im"], device=dev)[None],
                          torch.as_tensor(lab["Z"], device=dev)[None].double())
        cands["frame"] = (U0, V.solve_z(U0, t2, mask, args.lam, zr))
        cands[f"frame+adam{args.adam}"] = in_frame_adam(U0, t2, mask, args.lam, zr, args.adam, nocc=no)
        if args.no_baselines:
            cands = {"label": cands["label"]}       # label kept: per-molecule differences are reported against it
        for tag, model, a, chem, real in models:
            Ts = args.T or sorted({1, a.get("T", 1)})
            with torch.no_grad():
                traj = model(x["U0"], x["t2"], x["mask"], args.lam, x["zref"], T=max(Ts), chem=x["chem"])
            for t in Ts:
                U = traj[t][0].to(torch.complex128)
                cands[f"{tag}_t{t}"] = (U, V.solve_z(U, t2, mask, args.lam, zr))
                if args.net_adam:
                    cands[f"{tag}_t{t}+adam{args.net_adam}"] = in_frame_adam(U, t2, mask, args.lam, zr, args.net_adam, nocc=no)
        resid = {}
        for k, (U, Z) in cands.items():
            resid[k] = float(D.rel_residual(Z.double(), U.to(torch.complex128), t2)[0])
            tasks.append(((name, k), name, U[0].detach().cpu().numpy().astype(np.complex128),
                          Z[0].detach().cpu().numpy().astype(np.float64), t1))
        meta[name] = {"resid": resid, "norb": n}
        print(f"  {name:22s} " + " ".join(f"{k}:{v:.3f}" for k, v in resid.items()), flush=True)
    print(f"{len(tasks)} energies for {len(names)} molecules ({time.time()-t0:.0f}s to build)", flush=True)
    keys = list(next(iter(meta.values()))["resid"].keys())
    if args.skip_energy or args.dump:
        for k in keys:
            print(f"  {k:28s} resid {np.median([meta[n]['resid'][k] for n in meta]):.3f}")
        if args.dump:
            pickle.dump({"tasks": tasks, "meta": meta}, open(args.dump, "wb"))
            print(f"-> {args.dump}")
        return
    return run_energies(args, tasks, meta, names, keys, t0)


def run_energies(args, tasks, meta, names, keys, t0):
    res = {}
    with get_context("spawn").Pool(args.n_procs) as pool:
        for (name, k), E, cf, dt in pool.imap_unordered(energy_job, tasks):
            res.setdefault(name, {})[k] = {"E": E, "corr_frac": cf, "t": dt}
            print(f"  {name:22s} {k:24s} corr% {100*cf:7.1f}  resid {meta[name]['resid'][k]:.3f}  ({dt:.0f}s)", flush=True)
    print("\n=== median over molecules: % CCSD correlation energy recovered (t2 residual) ===")
    summ = {}
    has_lab = all("label" in res[n] for n in names)
    lab_cf = np.array([100 * res[n]["label"]["corr_frac"] for n in names]) if has_lab else None
    for k in keys:
        cf = np.array([100 * res[n][k]["corr_frac"] for n in names])
        rs = np.array([meta[n]["resid"][k] for n in names])
        diff = cf - lab_cf if has_lab else np.full_like(cf, np.nan)
        summ[k] = {"median_corr_pct": float(np.median(cf)), "mean_corr_pct": float(cf.mean()), "min_corr_pct": float(cf.min()),
                   "median_resid": float(np.median(rs)), "median_diff_vs_label": float(np.median(diff)),
                   "wins_vs_label": int((diff > 0).sum()), "n": len(cf)}
        print(f"  {k:28s} median {np.median(cf):6.1f}  mean {cf.mean():6.1f}  min {cf.min():7.1f}  resid {np.median(rs):.3f}"
              f"  vs label: median {np.median(diff):+5.1f}, better on {(diff > 0).sum()}/{len(cf)}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump({"summary": summ, "per_molecule": res, "meta": meta, "args": vars(args)}, open(ROOT / args.out, "w"), indent=1)
    print(f"-> {args.out}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
