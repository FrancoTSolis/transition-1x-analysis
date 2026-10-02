#!/usr/bin/env python3
"""Experiment 8: is the compressed-DF optimizer a (smooth) function of t2?

Trial runs showed that even from a gauge-canonical init, restarting ffsim's
L-BFGS from a phase-rotated copy of the same init lands on a different U
(phase-aligned distance ~ random) and, without regularization, a different Z.
Before generating labels at scale, quantify what controls this:

  variants (per molecule, all from the canonical init unless noted):
    base            reg=0,      maxiter=M
    long            reg=0,      maxiter=8*M            (convergence)
    reg             reg=0.005,  maxiter=M              (paper's lambda)
    reg_long        reg=0.005,  maxiter=8*M
    prox_<mu>       reg=0.005 + mu*||kappa - kappa_init||^2 (proximal/trust)
  each variant is run
    (i)  from the canonical init            -> (U, Z)
    (ii) from a random-phase-rotated init   -> (U', Z')      [gauge sensitivity]
    (iii) from the canonical init of t2 perturbed by 1% relative noise
                                             -> (U'', Z'')   [continuity]
  reporting resid, ||Z||, displacement from init, and for (ii)/(iii):
    phase-aligned ||U-U'||/sqrt(2n), ||Z-Z'||/||Z||, ||t2hat - t2hat'||/||t2||

Run from ccsd_amplitudes/ (ffsim venv, 1 core per job):
    python3 -m gauge_study.exp8_optimizer_determinism --n-mols 8 --n-procs 40
"""
from __future__ import annotations

import argparse
import os
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np  # noqa: E402

from .compressed_canonical import (canonical_exact_init, phase_dist,  # noqa: E402
                                   reconstruct_t2, rel_residual, z_mask)


def compress_prox(t2, Z0, U0, *, connectivity, regularization, znorm_ref, maxiter, prox_mu=0.0):
    """ffsim objective (+ regularizer) + optional proximal term on kappa params."""
    import jax
    import jax.numpy as jnp
    import scipy.optimize
    from opt_einsum import contract
    from ffsim.linalg.util import df_tensors_from_params, df_tensors_from_params_jax, \
        df_tensors_to_params
    from .compressed_canonical import diag_coulomb_indices_for
    jax.config.update("jax_enable_x64", True)
    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    n_tensors = Z0.shape[0]
    dci = diag_coulomb_indices_for(connectivity, norb)
    x0 = df_tensors_to_params(Z0, U0, dci)
    n_k = n_tensors * norb ** 2
    x0_k = jnp.asarray(x0[:n_k])
    t2_j = jnp.asarray(t2)

    def fun(x):
        Z, U = df_tensors_from_params_jax(x, n_tensors, norb, dci)
        rec = 1j * contract("mpq,map,mip,mbq,mjq->ijab", Z, U, U.conj(), U, U.conj(),
                            optimize="greedy")[:nocc, :nocc, nocc:, nocc:]
        loss = 0.5 * jnp.sum(jnp.abs(rec - t2_j) ** 2)
        if regularization:
            loss += regularization * jnp.abs(jnp.sum(jnp.abs(Z) ** 2) - znorm_ref)
        if prox_mu:
            loss += prox_mu * jnp.sum((x[:n_k] - x0_k) ** 2)
        return loss

    vg = jax.jit(jax.value_and_grad(fun))
    res = scipy.optimize.minimize(lambda x: tuple(map(np.asarray, vg(jnp.asarray(x)))),
                                  x0, method="L-BFGS-B", jac=True,
                                  options=dict(maxiter=maxiter))
    Z, U = df_tensors_from_params(res.x, n_tensors, norb, dci)
    return np.asarray(Z, float), np.asarray(U), int(res.nit), bool(res.success)


VARIANTS = {
    "base":      dict(reg=0.0,   mult=1, mu=0.0),
    "long":      dict(reg=0.0,   mult=8, mu=0.0),
    "reg":       dict(reg=0.005, mult=1, mu=0.0),
    "reg_long":  dict(reg=0.005, mult=8, mu=0.0),
    "prox1e-4":  dict(reg=0.005, mult=1, mu=1e-4),
    "prox1e-3":  dict(reg=0.005, mult=1, mu=1e-3),
    "prox1e-2":  dict(reg=0.005, mult=1, mu=1e-2),
}


def job(args):
    name, data_dir, conn, variant, maxiter, seed = args
    v = VARIANTS[variant]
    d = np.load(Path(data_dir) / f"{name}.npz")
    t2 = d["t2"].astype(np.float64)
    nocc, _, nvirt, _ = t2.shape
    norb = nocc + nvirt
    rng = np.random.default_rng(seed)
    init = canonical_exact_init(t2)
    m = z_mask(conn, norb)[None]
    Z0 = init.Z * m
    kw = dict(connectivity=conn, regularization=v["reg"], znorm_ref=init.znorm_full,
              maxiter=maxiter * v["mult"], prox_mu=v["mu"])
    t0 = time.time()
    Z, U, nit, ok = compress_prox(t2, Z0, init.U, **kw)
    t_opt = time.time() - t0
    # (ii) gauge-rotated init
    ph = np.exp(2j * np.pi * rng.random((2, norb)))
    Zg, Ug, _, _ = compress_prox(t2, Z0, init.U * ph[:, None, :], **kw)
    # (iii) perturbed input
    t2p = t2 * (1 + 0.01 * rng.standard_normal(t2.shape))
    t2p = 0.5 * (t2p + t2p.transpose(1, 0, 3, 2))          # keep ij<->ab symmetry
    initp = canonical_exact_init(t2p)
    Zp, Up, _, _ = compress_prox(t2p, initp.Z * m, initp.U, **{**kw, "znorm_ref": initp.znorm_full})
    nt2 = np.linalg.norm(t2)
    th = reconstruct_t2(Z, U, nocc)

    def cmp(Zb, Ub, t2b):
        return dict(
            dU=float(np.sqrt(sum(phase_dist(U[k], Ub[k]) ** 2 for k in range(2))) / np.sqrt(2 * norb)),
            dU_raw=float(np.linalg.norm(U - Ub) / np.sqrt(2 * norb)),
            dZ=float(np.linalg.norm(Z - Zb) / max(np.linalg.norm(Z), 1e-12)),
            dt2hat=float(np.linalg.norm(th - reconstruct_t2(Zb, Ub, nocc)) / nt2))
    out = dict(name=name, variant=variant, resid_init=rel_residual(t2, Z0, init.U),
               resid=rel_residual(t2, Z, U), znorm=float(np.linalg.norm(Z)),
               znorm_full=float(np.sqrt(init.znorm_full)), nit=nit, ok=ok, t=t_opt,
               disp=float(np.sqrt(sum(phase_dist(init.U[k], U[k]) ** 2 for k in range(2))) / np.sqrt(2 * norb)),
               gauge=cmp(Zg, Ug, t2), perturb=cmp(Zp, Up, t2p),
               dt2_in=float(np.linalg.norm(t2p - t2) / nt2),
               dinit_perturb=float(np.sqrt(sum(phase_dist(init.U[k], initp.U[k]) ** 2 for k in range(2))) / np.sqrt(2 * norb)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--names-file", default="gauge_study/names_n29_16_13.txt")
    ap.add_argument("--n-mols", type=int, default=8)
    ap.add_argument("--connectivity", default="square")
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--n-procs", type=int, default=40)
    ap.add_argument("--out", default="gauge_study/exp8_results.npy")
    args = ap.parse_args()
    names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    rng = np.random.default_rng(0)
    names = [names[i] for i in rng.choice(len(names), args.n_mols, replace=False)]
    jobs = [(n, args.data_dir, args.connectivity, v, args.maxiter, 100 + i)
            for i, n in enumerate(names) for v in args.variants]
    print(f"{len(jobs)} jobs ({args.n_mols} molecules x {len(args.variants)} variants), "
          f"connectivity={args.connectivity}, maxiter={args.maxiter}", flush=True)
    t0 = time.time()
    results = []
    with Pool(args.n_procs) as pool:
        for r in pool.imap_unordered(job, jobs):
            results.append(r)
            print(f"  [{len(results)}/{len(jobs)}] {r['name']} {r['variant']:9s} resid {r['resid_init']:.3f}->{r['resid']:.3f}"
                  f" |Z| {r['znorm']:.2f} nit {r['nit']} disp {r['disp']:.2f} | gauge dU {r['gauge']['dU']:.3f}"
                  f" dZ {r['gauge']['dZ']:.3f} dt2hat {r['gauge']['dt2hat']:.3f} | perturb dU {r['perturb']['dU']:.3f}"
                  f" dZ {r['perturb']['dZ']:.3f} dt2hat {r['perturb']['dt2hat']:.3f} ({time.time()-t0:.0f}s)", flush=True)
    np.save(args.out, np.array(results, dtype=object), allow_pickle=True)

    print(f"\n=== summary over {args.n_mols} molecules (medians), connectivity={args.connectivity} ===")
    print(f"{'variant':10s} {'resid':>6s} {'|Z|':>6s} {'|Zfull|':>7s} {'nit':>5s} {'conv':>5s} {'disp':>5s} | "
          f"{'gauge dU':>8s} {'dZ':>6s} {'dt2hat':>7s} | {'perturb dU':>10s} {'dZ':>6s} {'dt2hat':>7s} {'(dt2 in)':>8s} {'dinit':>6s}")
    for v in args.variants:
        rs = [r for r in results if r["variant"] == v]
        med = lambda f: np.median([f(r) for r in rs])  # noqa: E731
        print(f"{v:10s} {med(lambda r: r['resid']):6.3f} {med(lambda r: r['znorm']):6.2f} "
              f"{med(lambda r: r['znorm_full']):7.2f} {med(lambda r: r['nit']):5.0f} "
              f"{np.mean([r['ok'] for r in rs]):5.2f} {med(lambda r: r['disp']):5.2f} | "
              f"{med(lambda r: r['gauge']['dU']):8.3f} {med(lambda r: r['gauge']['dZ']):6.3f} "
              f"{med(lambda r: r['gauge']['dt2hat']):7.3f} | "
              f"{med(lambda r: r['perturb']['dU']):10.3f} {med(lambda r: r['perturb']['dZ']):6.3f} "
              f"{med(lambda r: r['perturb']['dt2hat']):7.3f} {med(lambda r: r['dt2_in']):8.3f} "
              f"{med(lambda r: r['dinit_perturb']):6.3f}")
    print(f"  init residual (masked exact DF): {np.median([r['resid_init'] for r in results]):.3f}")
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
