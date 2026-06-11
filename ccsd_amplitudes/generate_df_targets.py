#!/usr/bin/env python3
"""Generate double-factorization (U, Z) targets from cached RHF t2.

For each rhf_dataset/<name>.npz, runs ffsim's UCJOpSpinBalanced.from_t_amplitudes
(unconstrained, n_reps) and writes rhf_targets/<name>.npz with:
    Z          (n_reps, norb, norb)  diagonal-Coulomb matrices
    kappa_real (n_reps, norb, norb)  anti-symmetric  (= Re log U)
    kappa_imag (n_reps, norb, norb)  symmetric       (= Im log U)

These are used ONLY as a warm-start / anti-collapse supervision signal for the
reconstruction model: pure reconstruction from random init falls into a Z=0
local minimum (loss=1.0). Supervising Z (which is directly learnable, R^2~0.79)
for the first epochs moves the model off Z=0; reconstruction then refines U
gauge-free. `optimize=False` (exact truncated DF) is enough and fast.

Shardable:  --shard i --n-shards N
"""
import argparse
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--out-dir", default="rhf_targets")
    ap.add_argument("--n-reps", type=int, default=2)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    args = ap.parse_args()

    import ffsim
    from scipy.linalg import logm

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    files = sorted(data_dir.glob("*.npz"))
    files = [f for f in files if f.stem != "_index"]
    files = files[args.shard::args.n_shards]

    ok = skip = err = 0
    for i, f in enumerate(files):
        out_path = out_dir / f"{f.stem}.npz"
        if out_path.exists():
            skip += 1
            continue
        try:
            d = np.load(f)
            t2 = d["t2"].astype(float)
            t1 = d["t1"].astype(float)
            op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
                t2=t2, t1=t1, n_reps=args.n_reps, optimize=False,
            )
            U = op.orbital_rotations            # (n_reps, norb, norb) complex
            Z = op.diag_coulomb_mats[:, 0]       # (n_reps, norb, norb) real
            kr = np.stack([logm(U[k]).real for k in range(U.shape[0])])
            ki = np.stack([logm(U[k]).imag for k in range(U.shape[0])])
            np.savez(
                out_path,
                Z=Z.astype(np.float32),
                kappa_real=kr.astype(np.float32),
                kappa_imag=ki.astype(np.float32),
            )
            ok += 1
        except Exception as e:
            err += 1
            if err <= 5:
                print(f"  ERR {f.stem}: {type(e).__name__}: {e}", flush=True)
        if (i + 1) % 200 == 0:
            print(f"[shard {args.shard}] {i+1}/{len(files)}  ok={ok} skip={skip} err={err}",
                  flush=True)

    print(f"[shard {args.shard}] DONE ok={ok} skip={skip} err={err}", flush=True)


if __name__ == "__main__":
    main()
