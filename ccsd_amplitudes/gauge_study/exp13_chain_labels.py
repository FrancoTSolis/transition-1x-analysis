#!/usr/bin/env python3
"""Experiment 13: continuation ("chained") labels -- can the basin choice be made
consistent across neighbouring molecules without losing quality?

If the compressed-DF minima form a connected valley, warm-starting molecule k's
optimization from molecule k-1's *solution* (nearest neighbour in t2) should
keep the label on the same branch, giving a locally smooth label field at the
optimum's quality -- i.e. a lossless *and* regressable target.  Test: greedy
nearest-neighbour chains through the n29 group; compare chained labels with
the independent canonical labels (rhf_targets_compressed/square_reg0.005) on
  residual (quality) and consecutive-molecule dU (phase-aligned), dZ, dt2hat
  (smoothness), against the input dt2.

Run from ccsd_amplitudes/ (ffsim venv):
    python3 -m gauge_study.exp13_chain_labels --chains 4 --length 8 --n-procs 4
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

from .compressed_canonical import (canonical_exact_init, compress_from_init, phase_dist,  # noqa: E402
                                   reconstruct_t2, rel_residual, z_mask)


def _pin(cores):
    import multiprocessing
    if cores:
        wid = multiprocessing.current_process()._identity[0] - 1
        try:
            os.sched_setaffinity(0, {cores[wid % len(cores)]})
        except Exception:  # noqa: BLE001
            pass


def run_chain(args):
    names, data_dir, lab_dir, lam, maxiter = args
    out = []
    prev = None
    for k, name in enumerate(names):
        d = np.load(Path(data_dir) / f"{name}.npz")
        t2 = d["t2"].astype(np.float64)
        nocc, _, nvirt, _ = t2.shape
        norb = nocc + nvirt
        init = canonical_exact_init(t2)
        m = z_mask("square", norb)[None]
        if prev is None:
            Z0, U0 = init.Z * m, init.U
        else:
            Z0, U0 = prev
        t0 = time.time()
        r = compress_from_init(t2, Z0, U0, connectivity="square", regularization=lam,
                               znorm_ref=init.znorm_full, maxiter=maxiter)
        rec = dict(name=name, k=k, resid_init=rel_residual(t2, init.Z * m, init.U),
                   resid_chain=rel_residual(t2, r["Z"], r["U"]), nit=r["nit"], time=time.time() - t0,
                   U=r["U"], Z=r["Z"], t2hat=reconstruct_t2(r["Z"], r["U"], nocc), t2=t2)
        lab = Path(lab_dir) / f"{name}.npz"
        if lab.exists():
            c = np.load(lab)
            rec["U_fresh"] = c["U_re"] + 1j * c["U_im"]; rec["Z_fresh"] = c["Z"]
            rec["resid_fresh"] = float(c["resid"])
            rec["t2hat_fresh"] = reconstruct_t2(c["Z"], rec["U_fresh"], nocc)
        out.append(rec)
        prev = (r["Z"], r["U"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--lab-dir", default="rhf_targets_compressed/square_reg0.005")
    ap.add_argument("--names-file", default="gauge_study/names_n29_16_13.txt")
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--length", type=int, default=8)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--n-procs", type=int, default=4)
    ap.add_argument("--out", default="gauge_study/exp13_chain_labels.json")
    args = ap.parse_args()
    names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    names = [n for n in names if (Path(args.lab_dir) / f"{n}.npz").exists()]
    rng = np.random.default_rng(3)
    # greedy nearest-neighbour chains in raw t2 space
    X = np.stack([np.load(Path(args.data_dir) / f"{n}.npz")["t2"].astype(np.float32).ravel() for n in names])
    sq = (X ** 2).sum(1)
    chains = []
    used = set()
    for c in range(args.chains):
        cur = int(rng.integers(len(names)))
        while cur in used:
            cur = int(rng.integers(len(names)))
        chain = [cur]; used.add(cur)
        for _ in range(args.length - 1):
            d = sq[cur] + sq - 2 * X @ X[cur]
            d[list(used)] = np.inf
            cur = int(np.argmin(d)); chain.append(cur); used.add(cur)
        chains.append([names[i] for i in chain])
    print(f"{args.chains} chains x {args.length} molecules", flush=True)
    cores = sorted(os.sched_getaffinity(0))
    t0 = time.time()
    with Pool(args.n_procs, initializer=_pin, initargs=(cores,)) as pool:
        results = pool.map(run_chain, [(ch, args.data_dir, args.lab_dir, args.lam, args.maxiter) for ch in chains])
    print(f"optimizations done ({time.time()-t0:.0f}s)\n")

    rows = []
    print("=== per chain: consecutive-molecule distances (phase-aligned dU / ||U||, dt2hat / ||t2||) ===")
    for ci, res in enumerate(results):
        print(f"chain {ci}:")
        for k in range(len(res)):
            r = res[k]
            line = (f"  k={k} {r['name']:22s} resid init {r['resid_init']:.3f} | chained {r['resid_chain']:.3f}"
                    f" | fresh {r.get('resid_fresh', float('nan')):.3f} | nit {r['nit']}")
            if k > 0:
                p = res[k - 1]
                n = r["U"].shape[-1]
                dt2 = np.linalg.norm(r["t2"] - p["t2"]) / np.linalg.norm(r["t2"])
                dU_c = np.sqrt(sum(phase_dist(p["U"][j], r["U"][j]) ** 2 for j in range(2))) / np.sqrt(2 * n)
                dth_c = np.linalg.norm(r["t2hat"] - p["t2hat"]) / np.linalg.norm(r["t2"])
                dZ_c = np.linalg.norm(r["Z"] - p["Z"]) / max(np.linalg.norm(r["Z"]), 1e-12)
                line += f" || dt2 {dt2:.3f} | chained dU {dU_c:.3f} dZ {dZ_c:.3f} dt2hat {dth_c:.3f}"
                if "U_fresh" in r and "U_fresh" in p:
                    dU_f = np.sqrt(sum(phase_dist(p["U_fresh"][j], r["U_fresh"][j]) ** 2 for j in range(2))) / np.sqrt(2 * n)
                    dth_f = np.linalg.norm(r["t2hat_fresh"] - p["t2hat_fresh"]) / np.linalg.norm(r["t2"])
                    dZ_f = np.linalg.norm(r["Z_fresh"] - p["Z_fresh"]) / max(np.linalg.norm(r["Z_fresh"]), 1e-12)
                    line += f" | fresh dU {dU_f:.3f} dZ {dZ_f:.3f} dt2hat {dth_f:.3f}"
                    rows.append(dict(dt2=dt2, dU_c=dU_c, dZ_c=dZ_c, dth_c=dth_c, dU_f=dU_f, dZ_f=dZ_f, dth_f=dth_f,
                                     resid_c=r["resid_chain"], resid_f=r["resid_fresh"], resid_i=r["resid_init"]))
            print(line, flush=True)
    if rows:
        med = lambda k: np.median([x[k] for x in rows])  # noqa: E731
        print("\n=== medians over consecutive pairs ===")
        print(f"  input dt2/||t2||               {med('dt2'):.3f}")
        print(f"  residual: init {med('resid_i'):.3f}  chained {med('resid_c'):.3f}  fresh {med('resid_f'):.3f}")
        print(f"  dU  (phase-aligned, rel):  chained {med('dU_c'):.3f}   fresh {med('dU_f'):.3f}")
        print(f"  dZ  (rel):                 chained {med('dZ_c'):.3f}   fresh {med('dZ_f'):.3f}")
        print(f"  dt2hat (rel):              chained {med('dth_c'):.3f}   fresh {med('dth_f'):.3f}")
    json.dump(rows, open(args.out, "w"), indent=1)
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
