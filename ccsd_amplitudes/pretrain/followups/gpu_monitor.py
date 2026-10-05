#!/usr/bin/env python3
"""Log the real (nvidia-smi) GPU memory footprint of every compute process on one GPU, once per `--every` seconds.

One json line per sample: {"t", "gpu_used_MiB", "gpu_util", "procs": {pid: MiB}, "fts": {pid: MiB}} where "fts" are
the processes of the current user (owner from ps).  Read-only: it never touches any process.

Usage (on the GPU host):  python3 pretrain/followups/gpu_monitor.py --gpu 4 --out <log.jsonl> [--every 1] [--hours 30]
"""
from __future__ import annotations

import argparse
import getpass
import json
import subprocess
import time


def _q(args):
    return subprocess.run(["nvidia-smi", *args, "--format=csv,noheader,nounits"], capture_output=True, text=True,
                          timeout=30).stdout.strip().splitlines()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--every", type=float, default=1.0)
    ap.add_argument("--hours", type=float, default=30.0, help="stop after this many hours")
    args = ap.parse_args()
    me = getpass.getuser()
    uuid = None
    for ln in _q(["--query-gpu=index,uuid"]):
        i, u = [x.strip() for x in ln.split(",")]
        if int(i) == args.gpu:
            uuid = u
    assert uuid, f"GPU {args.gpu} not found"
    owners = {}
    t_end = time.time() + args.hours * 3600
    with open(args.out, "a") as f:
        while time.time() < t_end:
            t = time.time()
            try:
                g = _q([f"--id={args.gpu}", "--query-gpu=memory.used,utilization.gpu"])[0].split(",")
                procs = {}
                for ln in _q(["--query-compute-apps=pid,gpu_uuid,used_memory"]):
                    if not ln.strip():
                        continue
                    pid, gu, mem = [x.strip() for x in ln.split(",")]
                    if gu != uuid:
                        continue
                    procs[pid] = int(mem)
                    if pid not in owners:
                        owners[pid] = subprocess.run(["ps", "-o", "user=", "-p", pid], capture_output=True,
                                                     text=True).stdout.strip()
                rec = {"t": t, "gpu_used_MiB": int(g[0]), "gpu_util": int(g[1]), "procs": procs,
                       "fts": {p: m for p, m in procs.items() if owners.get(p) == me}}
                f.write(json.dumps(rec) + "\n")
                f.flush()
            except Exception as e:  # noqa: BLE001
                f.write(json.dumps({"t": t, "error": f"{type(e).__name__}: {e}"}) + "\n")
                f.flush()
            time.sleep(max(0.0, args.every - (time.time() - t)))


if __name__ == "__main__":
    main()
