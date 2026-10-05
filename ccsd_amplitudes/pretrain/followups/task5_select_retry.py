#!/usr/bin/env python3
"""Task 5: list the rl4n29f (molecule, chi) timing jobs whose fastest state build without concurrent own builds
was taken while the GPU was busy (another user's job) -> job list for a repeat pass (stdout, comma-separated)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
recs = [json.loads(ln) for ln in open(ROOT / "pretrain/opt_true/results/followups/task5_dm_scaling/bench.jsonl")]
recs = [r for r in recs if "error" not in r and r["cand"] == "rl4n29f"]
mon = [json.loads(ln) for ln in open(ROOT / "rl_runs/followups/task5_dm_scaling/gpu4_monitor_partA.jsonl")]
mon = [m for m in mon if "procs" in m]
ts = np.array([m["t"] for m in mon])
for r in recs:
    r["s1"] = r["t_start"] + r["t_state"]
for r in recs:
    r["overlap"] = sum(max(0.0, min(r["s1"], o["s1"]) - max(r["t_start"], o["t_start"])) for o in recs if o is not r) / r["t_state"]
    i0, i1 = np.searchsorted(ts, r["t_start"]), np.searchsorted(ts, r["s1"])
    u = np.array([m["gpu_util"] for m in mon[i0:i1]])
    r["busy"] = float(np.mean(u >= 95)) if len(u) else 1.0
best = {}
for r in recs:
    k = (r["name"], r["chi"])
    if r["overlap"] < 0.5 and r["busy"] <= float(sys.argv[1] if len(sys.argv) > 1 else 0.2):
        if k not in best or r["t_state"] < best[k]["t_state"]:
            best[k] = r
want = sorted({(r["name"], r["chi"], r["norb"]) for r in recs if r["chi"] <= 256}, key=lambda x: (x[2], x[0], x[1]))
todo = [(n, c) for n, c, _ in want if (n, c) not in best]
for n, c, nb in want:
    b = best.get((n, c))
    print(f"# n{nb} {n} chi{c}: " + (f"clean t_state {b['t_state']:.1f} (busy {b['busy']:.2f})" if b else "NO clean run"),
          file=sys.stderr)
print(",".join(f"{n}:rl4n29f:{c}" for n, c in todo))
