#!/usr/bin/env python3
"""Precompute the canonical exact-DF init (gauge_study.compressed_canonical.canonical_exact_init) for every
molecule in rhf_dataset/, grouped by shape, for multi-size optimize=True pretraining.

Output: rhf_inits_canonical/<nocc>_<nvirt>.npz with arrays
    names (M,), U_re/U_im (M,2,n,n) float32, Z (M,2,n,n) float32 (unmasked), znorm_full (M,), lam0 (M,), gap (M,)

Usage: python3 -m pretrain.opt_true.precompute_inits --n-procs 16
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from gauge_study.compressed_canonical import canonical_exact_init  # noqa: E402


def one(name):
    try:
        t2 = np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t2"].astype(np.float64)
        ini = canonical_exact_init(t2)
        return name, ini.U.astype(np.complex64), ini.Z.astype(np.float32), float(ini.znorm_full), float(ini.lam0), float(ini.gap)
    except Exception as e:  # noqa: BLE001
        return name, None, str(e), None, None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-procs", type=int, default=16)
    ap.add_argument("--out-dir", default="rhf_inits_canonical")
    args = ap.parse_args()
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    groups = defaultdict(list)
    for k, (n, no, nv) in idx.items():
        groups[(no, nv)].append(k)
    out = ROOT / args.out_dir
    out.mkdir(exist_ok=True)
    t0 = time.time()
    with Pool(args.n_procs) as pool:
        for (no, nv), names in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            f = out / f"{no}_{nv}.npz"
            if f.exists():
                continue
            res = pool.map(one, sorted(names), chunksize=8)
            ok = [r for r in res if r[1] is not None]
            bad = [r for r in res if r[1] is None]
            np.savez(f, names=np.array([r[0] for r in ok]), U_re=np.stack([r[1].real for r in ok]),
                     U_im=np.stack([r[1].imag for r in ok]), Z=np.stack([r[2] for r in ok]),
                     znorm_full=np.array([r[3] for r in ok]), lam0=np.array([r[4] for r in ok]),
                     gap=np.array([r[5] for r in ok]))
            print(f"  ({no},{nv}) n={no+nv}: {len(ok)} ok, {len(bad)} failed  ({time.time()-t0:.0f}s)", flush=True)
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
