#!/usr/bin/env python3
"""Score task-2 parameter sets with the Task-1 protocol (Lin et al.'s state-vector QSCI), identical to
pretrain/followups/task1_sqd_lin.py so the numbers can be put next to its candidates:

  exact LUCJ state (GPU, complex64 -> complex128, normalized); ffsim.sample_state_vector(1_000_000 shots,
  seed = default_rng(0)); uniformly random subset of 100_000 (default_rng(12345)); first iteration of
  diagonalize_fermionic_hamiltonian (10 batches x 4000, max_dim 4000, symmetrize_spin) via
  task2_common.task1_ci_strings (a copy of task1's lin_ci_strings, seed 0); qiskit_addon_sqd solve_sci with spin_sq = 0.0 per batch (identical batches solved once) in a
  torch-free process; mean / min / max over the 10 batches.  The lowest batch's eigenvector is re-checked by its
  full-space Rayleigh quotient on the GPU engine (d_full).

Parameter sets: the start points (label, rl4n29f) of every molecule that has a result, and x_best of every
result JSON under results/<objective>/ (optionally filtered).  Appends to results/score_task1_protocol.jsonl
(skips entries already scored).

  CUDA_VISIBLE_DEVICES=3 python -m pretrain.followups.task2_score [--only qsci] [--threads 3]
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
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--max-mem-gb", type=float, default=3.5)
    ap.add_argument("--only", nargs="*", default=None, help="objectives to include (lucj qsci); default both")
    ap.add_argument("--names", nargs="*", default=None)
    ap.add_argument("--no-starts", action="store_true")
    ap.add_argument("--starts-for", nargs="*", default=None, help="score the start points of these molecules too")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    out_path = C.RESULTS / "score_task1_protocol.jsonl"
    done = set()
    if out_path.exists():
        for ln in open(out_path):
            r = json.loads(ln)
            done.add((r["name"], r["cand"]))

    items = {}
    for obj in (args.only or ["lucj", "qsci"]):
        for f in sorted((C.RESULTS / obj).glob("*.json")):
            r = json.load(open(f))
            if args.names and r["name"] not in args.names:
                continue
            if "__tune" in f.stem or "smoke" in f.stem:
                continue
            z = np.load(f.with_suffix(".npz"))
            cand = f"{obj}:{r['start']}:{r['optimizer']}" + (f":{r['args']['tag']}" if r["args"].get("tag") else "")
            items.setdefault(r["name"], []).append((cand, z["x_best"]))
    for name in (args.starts_for or []):
        items.setdefault(name, [])
    if not args.no_starts:
        for name in list(items):
            for start in ("label", "rl4n29f"):
                U, Z = C.start_uz(name, start)
                items[name].insert(0, (f"start:{start}", C.uz_to_x(U, Z)))
    todo = [(n, c, x) for n, lst in items.items() for c, x in lst if (n, c) not in done]
    print(f"{len(todo)} parameter sets to score", flush=True)
    if args.dry_run:
        for n, c, _ in todo:
            print("  ", n, c)
        return
    by_mol = {}
    for n, c, x in todo:
        by_mol.setdefault(n, []).append((c, x))
    cct_all = json.load(open(C.ROOT / "pretrain/opt_true/results/ccsd_t_refs.json"))
    for name, lst in by_mol.items():
        ham = C.load_ham(name)
        cct = cct_all.get(name, {})
        log = C.LOGS / "score"
        log.mkdir(parents=True, exist_ok=True)
        rem = Remote(name, args.threads, args.max_mem_gb, 0, log / f"{name}.worker.log", sci=True)
        try:
            for cand, x in lst:
                s = rem.call({"cmd": "state_ci", "x": x})
                r = rem.call({"cmd": "solve", "ci": s["ci"], "spin_sq": C.T1_PROTOCOL["spin_sq"]})
                Es = np.array(r["E_batches"])
                dims = [int(len(sa)) for sa, _ in s["ci"]]
                eh, ec = ham["e_hf"], ham["e_ccsd"]
                rec = dict(name=name, cand=cand, norb=ham["norb"], e_hf=eh, e_ccsd=ec, e_ccsd_t=cct.get("e_ccsd_t"),
                           E_var=s["E_var"], corr_var=C.corr_pct(s["E_var"], eh, ec), entropy=s["entropy"],
                           p_hf=s["p_hf"], n_unique_1M=s["n_unique_1M"], n_unique_100k=s["n_unique_100k"], dims=dims,
                           n_distinct_batches=r["n_distinct"], E_sqd=Es.tolist(), E_sqd_mean=float(Es.mean()),
                           E_sqd_min=float(Es.min()), E_sqd_max=float(Es.max()),
                           corr_sqd_mean=C.corr_pct(Es.mean(), eh, ec), corr_sqd_min=C.corr_pct(Es.min(), eh, ec),
                           corr_sqd_max=C.corr_pct(Es.max(), eh, ec),
                           err_sqd_mean_vs_ccsdt_mHa=None if not cct else 1000 * (Es.mean() - cct["e_ccsd_t"]),
                           d_full=r.get("d_full"), asym=r.get("asym"), t_state_ci=s["t"], t_solve=r["t_sci"],
                           threads=args.threads, protocol=dict(C.T1_PROTOCOL),
                           host=os.uname().nodename, finished=time.strftime("%Y-%m-%d %H:%M:%S"))
                with open(out_path, "a") as fh:
                    fh.write(json.dumps(rec) + "\n")
                print(f"  {name:20s} {cand:30s} var {rec['corr_var']:6.2f}%  SQD mean {rec['corr_sqd_mean']:6.2f}% "
                      f"[{rec['corr_sqd_min']:.2f}, {rec['corr_sqd_max']:.2f}]  uniq100k {s['n_unique_100k']}  "
                      f"dims {min(dims)}-{max(dims)}  d_full {rec['d_full']:.1e}  ({s['t']:.0f}+{r['t_sci']:.0f}s)",
                      flush=True)
        finally:
            rem.close()


if __name__ == "__main__":
    main()
