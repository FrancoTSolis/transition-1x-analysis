#!/usr/bin/env python3
"""Tensor-network (MPS) LUCJ energies for pickled candidate tasks (policy_dump / eval_slot_energy format), CPU pool.

Writes the eval_slot_energy JSON layout (per_molecule[name][candidate] = {E, corr_frac}, meta), so
pretrain.opt_true.energy_table can summarise it.  Use a high bond dimension (256) for reported absolute energies.

Usage: python3 -m pretrain.rl.tn_eval --dump runs_ot/energy_tasks/n29val.pkl --chi 256 --n-workers 100 \
           --out pretrain/opt_true/results/energy_n29val_chi256.json
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from multiprocessing import get_context
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.rl.grpo import _reward_job, _worker_init  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="append", required=True)
    ap.add_argument("--chi", type=int, default=256)
    ap.add_argument("--n-workers", type=int, default=16)
    ap.add_argument("--basis-cache", default="rl_runs/tn_basis_cache")
    ap.add_argument("--stack-mem-gb", type=float, default=2.0)
    ap.add_argument("--skip-keys", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    tasks, meta = [], {}
    for f in args.dump:
        d = pickle.load(open(ROOT / f, "rb"))
        tasks += d["tasks"]
        for n, m in d["meta"].items():
            meta.setdefault(n, {"resid": {}, "norb": m["norb"]})["resid"].update(m["resid"])
    tasks = [t for t in {t[0]: t for t in tasks}.values() if t[0][1] not in args.skip_keys]
    kw = dict(max_bond=args.chi, basis_cache=str(ROOT / args.basis_cache), stack_mem=int(args.stack_mem_gb * (1 << 30)),
              block2_threads=1, tn_cache_items=1)
    jobs = [(k, n, U, Z, t1, "square", "tn", kw) for (k, n, U, Z, t1) in tasks]
    t0 = time.time()
    ham = {n: np.load(ROOT / "rhf_hamiltonians" / f"{n}.npz") for (_, n, *_r) in tasks}
    per = {}
    with get_context("spawn").Pool(args.n_workers, initializer=_worker_init, initargs=(str(ROOT / "rhf_hamiltonians"), 1)) as pool:
        for key, E, info in pool.imap_unordered(_reward_job, jobs):
            n, cand = key
            e_hf, e_cc = float(ham[n]["e_hf"]), float(ham[n]["e_ccsd"])
            per.setdefault(n, {})[cand] = {"E": E, "corr_frac": (e_hf - E) / (e_hf - e_cc), "t": info.get("t"),
                                           "discarded": info.get("discarded_sum"), "error": info.get("error")}
            print(f"  {n:22s} {cand:10s} corr% {100*per[n][cand]['corr_frac']:7.2f}  ({info.get('t', 0):.0f}s)", flush=True)
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"per_molecule": per, "meta": meta, "args": vars(args), "wall_s": time.time() - t0}, open(out, "w"), indent=1)
    print(f"-> {args.out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
