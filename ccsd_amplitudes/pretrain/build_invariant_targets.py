#!/usr/bin/env python3
"""Precompute gauge-free invariant targets (lam0, v0) for the retargeted head.

For each rhf_dataset/<name>.npz, eigendecompose the reshaped t2 matrix
    T[(i,a),(j,b)] = t2[i,j,a,b]
and store the top eigenpair in rhf_inv_targets/<name>.npz:
    lam0  float32   signed eigenvalue with largest |lambda|
    v0    (nocc, nvirt) float32   unit eigenvector as a matrix, sign-fixed
                                  so its largest-|entry| element is positive
    gap   float32   (|lam0| - |lam1|) / |lam0|   (Davis-Kahan stability)

This is the sufficient statistic of the n_reps=2 exact-DF (U, Z) labels
(see docs/eigenvector_learnability_study.md section 2): the full (U, Z) is
recovered at inference by the deterministic quadrature+eigh construction.
Pure numpy - no ffsim, no gauge left in the label pipeline.

Usage:  python3 -m pretrain.build_invariant_targets --n-procs 16
"""
import argparse
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np


def process_one(args):
    src, dst = args
    try:
        d = np.load(src)
        t2 = d["t2"].astype(np.float64)
        nocc, _, nvirt, _ = t2.shape
        T = t2.transpose(0, 2, 1, 3).reshape(nocc * nvirt, nocc * nvirt)
        eigs, vecs = np.linalg.eigh(T)
        order = np.argsort(np.abs(eigs))[::-1]
        lam0 = eigs[order[0]]
        lam1 = eigs[order[1]] if len(eigs) > 1 else 0.0
        v0 = vecs[:, order[0]]
        v0 = v0 * np.sign(v0[np.argmax(np.abs(v0))])
        gap = (abs(lam0) - abs(lam1)) / abs(lam0) if lam0 != 0 else 0.0
        np.savez(
            dst,
            lam0=np.float32(lam0),
            v0=v0.reshape(nocc, nvirt).astype(np.float32),
            gap=np.float32(gap),
        )
        return float(lam0)
    except Exception as e:
        print(f"ERR {src.stem}: {type(e).__name__}: {e}", flush=True)
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--out-dir", default="rhf_inv_targets")
    ap.add_argument("--n-procs", type=int, default=16)
    args = ap.parse_args()

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    tasks = []
    for f in sorted(Path(args.data_dir).glob("*.npz")):
        if f.stem == "_index":
            continue
        dst = out_dir / f.name
        if not dst.exists():
            tasks.append((f, dst))
    print(f"{len(tasks)} molecules to process -> {out_dir}/", flush=True)

    t0 = time.time()
    lams = []
    with Pool(args.n_procs) as pool:
        for i, lam in enumerate(pool.imap_unordered(process_one, tasks, chunksize=32)):
            if lam is not None:
                lams.append(lam)
            if (i + 1) % 2000 == 0:
                print(f"  {i+1}/{len(tasks)}  ({time.time()-t0:.0f}s)", flush=True)

    lams = np.array(lams)
    if len(lams):
        print(f"\nDONE {len(lams)} ok in {time.time()-t0:.0f}s")
        print(f"lam0 stats: mean={lams.mean():.4f} std={lams.std():.4f} "
              f"min={lams.min():.4f} max={lams.max():.4f} "
              f"frac_negative={(lams < 0).mean():.3f}")


if __name__ == "__main__":
    main()
