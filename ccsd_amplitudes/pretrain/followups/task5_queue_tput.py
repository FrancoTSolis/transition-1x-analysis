#!/usr/bin/env python3
"""Task 5: reward throughput of the dm TN engine through the real queue (the RL configuration), per norb.

Submits one batch of `--per-size` energy tasks per norb (rl4n29f parameters of the task-5 molecules, repeated, in
random order -- as in a GRPO step, consecutive tasks of a worker are mostly different molecules, so most tasks
rebuild their block2 MPO), waits for the batch, then the next size.  Records the batch wall time (-> energies per
GPU-hour with the workers serving the queue) and every task's worker time.

    python3 pretrain/followups/task5_queue_tput.py --queue-root rl_queue/task6_tn128 --sizes 44,37,33,29 \
        --per-size 12 --out pretrain/opt_true/results/followups/task5_dm_scaling/queue_tput.json
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl import reward_queue as RQ  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queue-root", required=True)
    ap.add_argument("--params", default="rl_runs/followups/task5_dm_scaling/params_rl4n29f.pkl")
    ap.add_argument("--sizes", default="44,37,33,29")
    ap.add_argument("--per-size", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    d = pickle.load(open(ROOT / args.params, "rb"))
    rng = np.random.default_rng(args.seed)
    out = {"args": vars(args), "batches": []}
    for n in [int(x) for x in args.sizes.split(",")]:
        tasks = []
        for (key, name, U, Z, t1) in d["tasks"]:
            h = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
            if int(h["norb"]) != n:
                continue
            tasks.append({"name": name, "U": U, "Z": Z, "t1": t1, "norb": n,
                          "nelec": (int(h["nelec_a"]), int(h["nelec_b"])), "kind": "tn"})
        sub = [tasks[i % len(tasks)] for i in rng.permutation(args.per_size)]
        t0 = time.time()
        ids = RQ.submit(ROOT / args.queue_root, sub)
        res = RQ.collect(ROOT / args.queue_root, ids, stale_s=1800.0)
        wall = time.time() - t0
        per = [{"name": s["name"], "E": res[i][0], **res[i][1]} for s, i in zip(sub, ids)]
        rec = {"norb": n, "n_tasks": len(sub), "wall_s": wall, "energies_per_h": 3600 * len(sub) / wall,
               "s_per_energy": wall / len(sub), "task_t_median": float(np.median([p.get("t", np.nan) for p in per])),
               "n_failed": sum(1 for p in per if not np.isfinite(p["E"])), "tasks": per, "t_start": t0}
        out["batches"].append(rec)
        print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in rec.items() if k != "tasks"}),
              flush=True)
        Path(ROOT / args.out).parent.mkdir(parents=True, exist_ok=True)
        json.dump(out, open(ROOT / args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
