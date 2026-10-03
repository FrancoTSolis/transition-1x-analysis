#!/usr/bin/env python3
"""Energy tasks (for pretrain.opt_true.eval_slot_energy --from-dump) from GRPO slot policies on any molecule set.

Candidates per molecule: the deterministic mean action of each given policy checkpoint (tag), optionally the
optimize=True label from --labels-dir.  Run on a GPU or CPU; energies are computed elsewhere.

Usage: python3 -m pretrain.rl.policy_dump --names-file gauge_study/names_norb17.txt \
          --policy rl_runs/grpo_slot_v2/policy_best.pt --tag rl1 [--policy ... --tag ...] \
          --labels-dir rhf_targets_compressed_n17 --out runs_ot/energy_tasks/n17.pkl
"""
from __future__ import annotations

import argparse
import copy
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.opt_true import dftorch as D  # noqa: E402
from pretrain.opt_true.eval_slot import load_slot  # noqa: E402
from pretrain.rl.grpo_slot import Mol, Policy, to_flat  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names-file", required=True)
    ap.add_argument("--policy", action="append", default=[], help="GRPO policy checkpoint (policy_*.pt) or 'init:<slot ckpt>'")
    ap.add_argument("--tag", action="append", default=[])
    ap.add_argument("--labels-dir", default=None)
    ap.add_argument("--lam", type=float, default=0.005)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    dev = args.device
    names = [ln.strip() for ln in open(ROOT / args.names_file) if ln.strip()]
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    tasks, meta = [], {n: {"resid": {}, "norb": idx[n][0]} for n in names}
    for pol_path, tag in zip(args.policy, args.tag):
        if pol_path.startswith("init:"):
            init, prefix_T, sd = pol_path[5:], 0, None
        else:
            ck = torch.load(ROOT / pol_path, map_location=dev, weights_only=False)
            init, prefix_T, sd = ck["args"]["init"], ck["args"].get("prefix_T", 0), ck["policy"]
        slot, a_, _, _ = load_slot(ROOT / init, dev)
        frozen = copy.deepcopy(slot).eval() if prefix_T else None
        pol = Policy(slot, a_["d"]).to(dev)
        if sd is not None:
            pol.load_state_dict(sd)
        pol.eval()
        for n in names:
            m = Mol(n, idx, dev, args.lam, frozen, prefix_T)
            with torch.no_grad():
                dK, dZ = pol(m.x, args.lam)
                U, Z, r = m.realize(to_flat(dK.double(), dZ.double(), m.fi))
            tasks.append(((n, tag), n, U, Z, m.t1))
            meta[n]["resid"][tag] = r
        print(f"{tag}: median residual {np.median([meta[n]['resid'][tag] for n in names]):.3f}", flush=True)
    if args.labels_dir:
        for n in names:
            lab = np.load(ROOT / args.labels_dir / "square_reg0.005" / f"{n}.npz")
            U, Z = lab["U_re"] + 1j * lab["U_im"], lab["Z"]
            t1 = np.load(ROOT / "rhf_dataset" / f"{n}.npz")["t1"].astype(np.float64)
            tasks.append(((n, "label"), n, U, Z, t1))
            meta[n]["resid"]["label"] = float(lab["resid"])
    pickle.dump({"tasks": tasks, "meta": meta}, open(ROOT / args.out, "wb"))
    print(f"{len(tasks)} tasks -> {args.out}")


if __name__ == "__main__":
    main()
