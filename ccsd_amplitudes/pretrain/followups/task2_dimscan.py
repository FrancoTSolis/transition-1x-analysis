#!/usr/bin/env python3
"""Subspace energy at fixed dimension: the k most probable single-spin strings of the exact LUCJ state (exact
marginals p(a) = sum_b |psi[a, b]|^2 on the GPU), the k x k Cartesian-product subspace (as QSCI with
symmetrize_spin), diagonalized with qiskit_addon_sqd solve_sci (spin_sq = 0) in a torch-free process and checked by
the full-space Rayleigh quotient.  Separates "better configurations" from "more configurations": QSCI-optimized
states win the fixed-sample-count objective partly by spreading the distribution (larger subspaces).

  CUDA_VISIBLE_DEVICES=3 python -m pretrain.followups.task2_dimscan --names C2H3N_rxn2858_P --ks 150 300 450 600 800
Appends to results/dimscan.jsonl; parameter sets as in task2_score (start points + x_best of every result).
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.followups import task2_common as C  # noqa: E402
from pretrain.followups.task2_opt import Remote  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--ks", type=int, nargs="+", default=[150, 300, 450, 600, 800])
    ap.add_argument("--only", nargs="*", default=None, help="objectives of the x_best sets (lucj qsci)")
    ap.add_argument("--cands", nargs="*", default=None, help="restrict to these candidate labels")
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--max-mem-gb", type=float, default=2.0)
    args = ap.parse_args()
    out_path = C.RESULTS / "dimscan.jsonl"
    done = set()
    if out_path.exists():
        for ln in open(out_path):
            r = json.loads(ln)
            done.add((r["name"], r["cand"], r["k"]))
    for name in args.names:
        items = []
        for start in ("label", "rl4n29f"):
            U, Z = C.start_uz(name, start)
            items.append((f"start:{start}", C.uz_to_x(U, Z)))
        for obj in (args.only or ["lucj", "qsci"]):
            for f in sorted((C.RESULTS / obj).glob(f"{name}__*.json")):
                r = json.load(open(f))
                tag = r["args"].get("tag") or ""
                if tag.startswith(("tune", "smoke", "long")):
                    continue
                cand = f"{obj}:{r['start']}:{r['optimizer']}" + (f":{tag}" if tag else "")
                items.append((cand, np.load(f.with_suffix(".npz"))["x_best"]))
        if args.cands:
            items = [it for it in items if it[0] in args.cands]
        items = [it for it in items if any((name, it[0], k) not in done for k in args.ks)]
        if not items:
            continue
        ham = C.load_ham(name)
        fci = {}
        fp = C.ROOT / "pretrain/opt_true/results/followups/task1_sqd_vs_lin/fci_refs.json"
        if fp.exists():
            fci = json.load(open(fp)).get(name, {})
        log = C.LOGS / "dimscan"
        log.mkdir(parents=True, exist_ok=True)
        rem = Remote(name, args.threads, args.max_mem_gb, 0, log / f"{name}.worker.log", sci=True)
        try:
            for cand, x in items:
                m = rem.call({"cmd": "marginals", "x": x})
                order = np.argsort(-m["p_string"], kind="stable")
                for k in args.ks:
                    if (name, cand, k) in done:
                        continue
                    top = np.sort(m["strings"][order[:k]])
                    t0 = time.time()
                    r = rem.call({"cmd": "solve", "ci": [(top, top)], "spin_sq": C.T1_PROTOCOL["spin_sq"]})
                    E = r["E_batches"][0]
                    rec = dict(name=name, cand=cand, k=k, n_det=k * k, E=E, E_var=m["E_var"],
                               corr=C.corr_pct(E, ham["e_hf"], ham["e_ccsd"]),
                               corr_var=C.corr_pct(m["E_var"], ham["e_hf"], ham["e_ccsd"]),
                               err_fci_mHa=None if not fci else 1000 * (E - fci["e_fci"]),
                               p_covered=float(m["p_string"][order[:k]].sum()), d_full=r.get("d_full"),
                               t=time.time() - t0, finished=time.strftime("%Y-%m-%d %H:%M:%S"))
                    with open(out_path, "a") as fh:
                        fh.write(json.dumps(rec) + "\n")
                    print(f"  {name} {cand:24s} k {k:5d}  E {E:.8f}  corr {rec['corr']:6.2f}%  "
                          f"err_FCI {rec['err_fci_mHa'] if fci else float('nan'):7.2f} mHa  marg {rec['p_covered']:.4f}"
                          f"  d_full {rec['d_full']:.1e}  ({rec['t']:.0f}s)", flush=True)
        finally:
            rem.close()


if __name__ == "__main__":
    main()
