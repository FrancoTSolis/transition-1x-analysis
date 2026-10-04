#!/usr/bin/env python3
"""Exact energies for pickled candidate tasks (eval_slot_energy / policy_dump format) through the GPU reward queue.

Writes the same JSON layout as pretrain.opt_true.eval_slot_energy (per_molecule[name][candidate] = {E, corr_frac},
meta[name][resid]), so pretrain.opt_true.energy_table can summarise it.

Usage: python3 -m pretrain.rl.queue_eval --dump runs_ot/energy_tasks/largeval_base.pkl --queue-root rl_queue/main \
           --out pretrain/opt_true/results/energy_largeval_base.json
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.rl import reward_queue as RQ  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="append", required=True)
    ap.add_argument("--queue-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--skip-keys", nargs="*", default=[])
    ap.add_argument("--only-names", nargs="*", default=None, help="restrict to these molecules")
    ap.add_argument("--kind", default="exact", choices=["exact", "tn"], help="worker kind serving this queue root")
    args = ap.parse_args()
    tasks, meta = [], {}
    for f in args.dump:
        d = pickle.load(open(ROOT / f, "rb"))
        tasks += d["tasks"]
        for n, m in d["meta"].items():
            meta.setdefault(n, {"resid": {}, "norb": m["norb"]})["resid"].update(m["resid"])
    tasks = [t for t in {t[0]: t for t in tasks}.values() if t[0][1] not in args.skip_keys
             and (args.only_names is None or t[1] in args.only_names)]
    ham = {}
    for (_, n, *_r) in tasks:
        if n not in ham:
            h = np.load(ROOT / "rhf_hamiltonians" / f"{n}.npz")
            ham[n] = (int(h["norb"]), (int(h["nelec_a"]), int(h["nelec_b"])), float(h["e_hf"]), float(h["e_ccsd"]))
    sub = [{"name": n, "U": U, "Z": Z, "t1": t1, "norb": ham[n][0], "nelec": ham[n][1], "kind": args.kind}
           for (k, n, U, Z, t1) in tasks]
    t0 = time.time()
    ids = RQ.submit(ROOT / args.queue_root, sub)
    print(f"submitted {len(ids)} tasks", flush=True)
    # TN energies can run 30 min: be slow to declare a worker dead (its heartbeat comes from a child process)
    res = RQ.collect(ROOT / args.queue_root, ids, verbose=True, stale_s=1800.0 if args.kind == "tn" else 300.0)
    per = {}
    for (k, n, *_r), i in zip(tasks, ids):
        E, info = res[i]
        e_hf, e_cc = ham[n][2], ham[n][3]
        per.setdefault(n, {})[k[1]] = {"E": E, "corr_frac": (e_hf - E) / (e_hf - e_cc), "t": info.get("t"),
                                       "worker": info.get("worker"), "error": info.get("error")}
    nbad = sum(1 for n in per for k in per[n] if per[n][k]["error"])
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"per_molecule": per, "meta": meta, "args": vars(args), "wall_s": time.time() - t0}, open(out, "w"), indent=1)
    print(f"-> {args.out}  ({time.time()-t0:.0f}s, {nbad} errors)")


if __name__ == "__main__":
    main()
