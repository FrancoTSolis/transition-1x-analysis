#!/usr/bin/env python3
"""Wall time of one rl4n29f inference (the amortized alternative to per-molecule optimization).

Per molecule: (a) input build = Mol(...) of pretrain.rl.grpo_slot (CCSD-amplitude features from disk, frame,
3 frozen pretrained recycles), (b) final policy recycle + realize (U from the generator, Z by least squares) --
exactly what pretrain.rl.policy_dump does.  Median of --reps repeats after one warm-up.

  CUDA_VISIBLE_DEVICES=3 python -m pretrain.followups.task2_infer_time --device cuda
  OMP_NUM_THREADS=1 python -m pretrain.followups.task2_infer_time --device cpu --threads 1
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.followups import task2_common as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--policy", default="rl_runs/grpo_n29_tn/policy_last.pt")
    ap.add_argument("--names-file", default="pretrain/rl/small_val.txt")
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    from pretrain.opt_true.eval_slot import load_slot
    from pretrain.rl.grpo_slot import Mol, Policy, to_flat
    dev = args.device
    sync = (lambda: torch.cuda.synchronize()) if dev.startswith("cuda") else (lambda: None)
    t0 = time.time()
    ck = torch.load(C.ROOT / args.policy, map_location=dev, weights_only=False)
    init, prefix_T, sd = ck["args"]["init"], ck["args"].get("prefix_T", 0), ck["policy"]
    slot, a_, _, _ = load_slot(C.ROOT / init, dev)
    frozen = copy.deepcopy(slot).eval() if prefix_T else None
    pol = Policy(slot, a_["d"]).to(dev)
    pol.load_state_dict(sd)
    pol.eval()
    t_load = time.time() - t0
    idx = json.load(open(C.ROOT / "rhf_dataset" / "_index.json"))
    names = [ln.strip() for ln in open(C.ROOT / args.names_file) if ln.strip()]
    lam = 0.005
    rows = {}
    for n in names:
        tb, tp = [], []
        for rep in range(args.reps + 1):
            sync()
            t0 = time.perf_counter()
            with torch.no_grad():
                m = Mol(n, idx, dev, lam, frozen, prefix_T)
            sync()
            t1 = time.perf_counter()
            with torch.no_grad():
                dK, dZ = pol(m.x, lam)
                U, Z, r = m.realize(to_flat(dK.double(), dZ.double(), m.fi))
            sync()
            t2 = time.perf_counter()
            if rep:
                tb.append(t1 - t0)
                tp.append(t2 - t1)
        U0, Z0 = C.start_uz(n, "rl4n29f")
        rows[n] = dict(norb=idx[n][0], t_build_s=float(np.median(tb)), t_policy_s=float(np.median(tp)),
                       t_total_s=float(np.median(np.array(tb) + np.array(tp))),
                       max_dU_vs_dump=float(np.abs(U - U0).max()), max_dZ_vs_dump=float(np.abs(Z - Z0).max()))
        print(n, json.dumps(rows[n]), flush=True)
    dev_name = torch.cuda.get_device_name() if dev.startswith("cuda") else f"cpu x{args.threads} threads"
    out = dict(device=dev_name, t_load_model_s=t_load, reps=args.reps, per_molecule=rows,
               median_total_s=float(np.median([r["t_total_s"] for r in rows.values()])))
    C.dumpj(out, C.RESULTS / f"infer_time_{'gpu' if dev.startswith('cuda') else 'cpu'}.json")
    print(json.dumps({k: v for k, v in out.items() if k != "per_molecule"}))


if __name__ == "__main__":
    main()
