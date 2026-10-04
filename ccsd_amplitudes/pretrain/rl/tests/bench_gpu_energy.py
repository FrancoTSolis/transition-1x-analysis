#!/usr/bin/env python3
"""Benchmark pretrain/rl/gpu_energy.py: wall time per energy (with phase breakdown), peak GPU memory,
and the batched API (energies() on n_batch parameter sets of one molecule).

  CUDA_VISIBLE_DEVICES=2 python -m pretrain.rl.tests.bench_gpu_energy --names C3H4_rxn2391_P \
      --tasks runs_ot/energy_tasks/baselines.pkl --repeat 3 --batch 8 --out pretrain/rl/tests/results/bench_n16.json
Parameters come from a task pickle (same molecule) when available, else random (U Haar, Z ~ N(0, 0.3)).
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402

from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac  # noqa: E402


def haar(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def params_for(name, norb, k, tasks_file, n, rng):
    out = []
    if tasks_file:
        for f in tasks_file:
            for (nm, cand), _, U, Z, t1 in pickle.load(open(ROOT / f, "rb"))["tasks"]:
                if nm == name:
                    out.append((U, Z, t1))
    while len(out) < n:  # random parameters near the identity (like trained models)
        U = np.stack([haar(norb, rng) for _ in range(2)])
        Z = 0.3 * rng.normal(size=(2, norb, norb))
        Z = 0.5 * (Z + Z.transpose(0, 2, 1))
        out.append((U, Z, 0.05 * rng.normal(size=(k, norb - k))))
    return out[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--tasks", nargs="*", default=[])
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--dtype", default="c8")
    ap.add_argument("--fuse", type=int, nargs="+", default=[1])
    ap.add_argument("--max-mem-gb", type=float, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rng = np.random.default_rng(0)
    dt = torch.complex64 if args.dtype == "c8" else torch.complex128
    res = []
    for name in args.names:
        for fuse in args.fuse:
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
            t0 = time.time()
            eng = LUCJEnergyGPU.from_npz(name, dtype=dt, fuse=bool(fuse), max_mem_gb=args.max_mem_gb)
            torch.cuda.synchronize()
            t_setup = time.time() - t0
            P = params_for(name, eng.norb, eng.k, args.tasks, max(args.batch, args.repeat), rng)
            eng.energy(*P[0][:2], t1=P[0][2])                       # warm-up (numba JIT, caches)
            per = []
            for i in range(args.repeat):
                U, Z, t1 = P[i]
                E = eng.energy(U, Z, t1=t1)
                per.append(dict(eng.timing, E=E, cf=corr_frac(E, eng.e_hf, eng.e_ccsd)))
            torch.cuda.synchronize()
            t0 = time.time()
            Es = eng.energies([(U, Z, t1) for U, Z, t1 in P[: args.batch]])
            torch.cuda.synchronize()
            t_batch = time.time() - t0
            r = dict(name=name, norb=eng.norb, k=eng.k, dim=eng.dim, fuse=fuse, t_top=eng.t_top, dtype=args.dtype,
                     nq=eng.nq, K=eng.K, tile=eng.T, mem_plan=eng.mem_plan,
                     peak_gb=torch.cuda.max_memory_allocated() / 1e9,
                     peak_reserved_gb=torch.cuda.max_memory_reserved() / 1e9, t_setup=t_setup,
                     per_energy=per, t_energy_mean=float(np.mean([p["total"] for p in per])),
                     t_batch=t_batch, n_batch=args.batch, batch_E=Es,
                     gpu=torch.cuda.get_device_name())
            res.append(r)
            ph = {k: float(np.mean([p[k] for p in per])) for k in per[0] if k not in ("E", "cf")}
            print(f"{name} n={eng.norb} k={eng.k} dim={eng.dim} fuse={fuse} t={eng.t_top}: "
                  f"{r['t_energy_mean']:.2f} s/energy  "
                  + " ".join(f"{k}={v:.3f}" for k, v in ph.items())
                  + f"  | batch {args.batch}: {t_batch:.2f} s  | peak {r['peak_gb']:.2f} GB "
                  f"(reserved {r['peak_reserved_gb']:.2f})  setup {t_setup:.1f}s  tile {eng.T}", flush=True)
            eng.release()
            del eng
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
