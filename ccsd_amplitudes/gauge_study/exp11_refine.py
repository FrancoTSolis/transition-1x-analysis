#!/usr/bin/env python3
"""Experiment 11: value of the amortized init as a warm start (predict-then-refine).

For val molecules of the n29 group, run ffsim's compressed optimization
(square, lambda=0.005) from
  (a) the canonical exact-DF init,
  (b) the network's output (compressed_recon checkpoint, best-of-K heads),
recording the t2 residual after every L-BFGS iteration (callback).  Reports the
iteration count needed from each start to reach given residual levels and the
final residuals -- i.e. how much per-molecule optimization the network saves,
and whether refining the network output reaches the same quality as the label.

Run from ccsd_amplitudes/ with the TRAIN venv (torch + ffsim):
    python3 -m gauge_study.exp11_refine --ckpt checkpoints_comp_recon_n29/best.pt --n-mols 8
"""
from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402

from .compressed_canonical import canonical_exact_init, rel_residual, z_mask  # noqa: E402


def _pin(cores):
    import multiprocessing
    if cores:
        wid = multiprocessing.current_process()._identity[0] - 1
        try:
            os.sched_setaffinity(0, {cores[wid % len(cores)]})
        except Exception:  # noqa: BLE001
            pass


def refine_job(args):
    """Run L-BFGS from (Z0, U0) and record residual per iteration."""
    name, t2, Z0, U0, lam, znorm_ref, maxiter, tag = args
    import jax
    import jax.numpy as jnp
    import scipy.optimize
    from opt_einsum import contract
    from ffsim.linalg.util import df_tensors_from_params, df_tensors_from_params_jax, df_tensors_to_params
    from .compressed_canonical import diag_coulomb_indices_for
    jax.config.update("jax_enable_x64", True)
    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    dci = diag_coulomb_indices_for("square", norb)
    t2_j = jnp.asarray(t2)

    def fun(x):
        Z, U = df_tensors_from_params_jax(x, 2, norb, dci)
        rec = 1j * contract("mpq,map,mip,mbq,mjq->ijab", Z, U, U.conj(), U, U.conj(),
                            optimize="greedy")[:nocc, :nocc, nocc:, nocc:]
        loss = 0.5 * jnp.sum(jnp.abs(rec - t2_j) ** 2)
        loss += lam * jnp.abs(jnp.sum(jnp.abs(Z) ** 2) - znorm_ref)
        return loss

    vg = jax.jit(jax.value_and_grad(fun))
    hist = []

    def cb(xk):
        Z, U = df_tensors_from_params(xk, 2, norb, dci)
        hist.append(rel_residual(t2, np.asarray(Z, float), np.asarray(U)))

    x0 = df_tensors_to_params(Z0, U0, dci)
    t0 = time.time()
    res = scipy.optimize.minimize(lambda x: tuple(map(np.asarray, vg(jnp.asarray(x)))), x0,
                                  method="L-BFGS-B", jac=True, callback=cb, options=dict(maxiter=maxiter))
    return dict(name=name, tag=tag, resid0=rel_residual(t2, Z0, U0), hist=hist, nit=int(res.nit),
                time=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default="checkpoints_comp_recon_n29/best.pt")
    ap.add_argument("--hyps", type=int, default=1)
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--names-file", default="gauge_study/names_n29_16_13.txt")
    ap.add_argument("--n-mols", type=int, default=8)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--n-procs", type=int, default=16)
    ap.add_argument("--out", default="gauge_study/exp11_refine.json")
    args = ap.parse_args()

    import torch
    from pretrain.dataset import CCSDAmplitudeDataset
    from pretrain.model import ModelConfig, PretrainingModel
    from scipy.linalg import expm

    # same split as train.py (seed 42, val_split 0.1, names filter)
    ds = CCSDAmplitudeDataset(args.data_dir, with_init=True, connectivity="square")
    allowed = {ln.strip() for ln in open(args.names_file) if ln.strip()}
    gen = torch.Generator().manual_seed(42)
    idx = torch.randperm(len(ds), generator=gen).tolist()
    idx = [i for i in idx if ds.names[i] in allowed]
    n_val = int(len(idx) * 0.1)
    val_idx = idx[len(idx) - n_val:][:args.n_mols]

    cfg = ModelConfig(embed_dim=192, num_layers=6, num_heads=8, n_reps=2, dropout=0.0,
                      predict_residual=True, residual_hyps=args.hyps)
    model = PretrainingModel(cfg)
    sd = torch.load(args.ckpt, map_location="cpu", weights_only=False)["model_state_dict"]
    model.load_state_dict(sd, strict=False)
    model.eval()

    jobs = []
    with torch.no_grad():
        for i in val_idx:
            s = ds[i]
            b = CCSDAmplitudeDataset.collate_fn([s])
            out = model(b)
            no, nv = s["nocc"], s["nvirt"]; n = no + nv
            t2 = s["t2"].numpy().astype(np.float64)
            init = canonical_exact_init(t2)
            m = z_mask("square", n)[None]
            Z0 = init.Z * m
            hyps = out.get("residual_hyps") or [out]
            cands = []
            for h in hyps:
                kr = h["dkappa_re"][0].numpy().astype(np.float64)[:, :n, :n]
                ki = h["dkappa_im"][0].numpy().astype(np.float64)[:, :n, :n]
                dz = h["dz"][0].numpy().astype(np.float64)[:, :n, :n] * m[0]
                U = np.stack([init.U[k] @ expm(kr[k] + 1j * ki[k]) for k in range(2)])
                Z = Z0 + dz
                cands.append((rel_residual(t2, Z, U), Z, U))
            r_best, Zn, Un = min(cands, key=lambda c: c[0])
            jobs.append((s["name"], t2, Z0, init.U, args.lam, init.znorm_full, args.maxiter, "exact-DF init"))
            jobs.append((s["name"], t2, Zn, Un, args.lam, init.znorm_full, args.maxiter, "network init"))
    print(f"{len(val_idx)} molecules x 2 starts", flush=True)
    cores = sorted(os.sched_getaffinity(0))
    t0 = time.time()
    res = []
    with Pool(args.n_procs, initializer=_pin, initargs=(cores,)) as pool:
        for r in pool.imap_unordered(refine_job, jobs):
            res.append(r)
            print(f"  {r['name']:24s} {r['tag']:14s} resid0={r['resid0']:.3f} -> "
                  f"{[round(r['hist'][k-1], 3) if len(r['hist']) >= k else None for k in (10, 25, 50, 100, 200)]} final {r['hist'][-1]:.3f} "
                  f"(nit {r['nit']}, {r['time']:.0f}s)", flush=True)
    json.dump(res, open(args.out, "w"))
    print("\n=== median residual vs L-BFGS iterations (square, lambda=%.3f) ===" % args.lam)
    ks = [0, 10, 25, 50, 100, 200, 500]
    for tag in ("exact-DF init", "network init"):
        rs = [r for r in res if r["tag"] == tag]
        row = []
        for k in ks:
            vals = [r["resid0"] if k == 0 else (r["hist"][min(k, len(r["hist"])) - 1] if r["hist"] else r["resid0"]) for r in rs]
            row.append(f"{np.median(vals):.3f}")
        print(f"  {tag:14s} " + "  ".join(f"it{k}={v}" for k, v in zip(ks, row)))
    for thr in (0.6, 0.5, 0.45, 0.4):
        for tag in ("exact-DF init", "network init"):
            rs = [r for r in res if r["tag"] == tag]
            its = [next((i + 1 for i, v in enumerate(r["hist"]) if v <= thr), None) for r in rs]
            reached = [i for i in its if i is not None]
            print(f"  iterations to reach resid<={thr}: {tag:14s} median {np.median(reached) if reached else 'n/a'}  "
                  f"(reached by {len(reached)}/{len(rs)})")
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
