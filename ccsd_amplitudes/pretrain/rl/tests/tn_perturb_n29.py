#!/usr/bin/env python3
"""GRPO-like perturbation group around n29 parameter sets (no exact energies exist at norb 29).
Same perturbation model as tn_perturbed_refs.py; writes <out>.npz/.json in the same format (E = null), so that
tn_bench_split.py --perturbed can evaluate it at several bond dimensions (self-consistency of the ranking).
Usage: python3 pretrain/rl/tests/tn_perturb_n29.py --params-dir <dir with bench_params_<name>.npz> \
          --names C3H5N3_rxn2003_P --tag net --group 8 --out pretrain/rl/tests/results/perturbed_n29
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pretrain.rl.tests.tn_perturbed_refs import perturb, polar  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params-dir", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--tag", default="net")
    ap.add_argument("--group", type=int, default=8)
    ap.add_argument("--seed", type=int, default=4321)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    names = args.names.split(",")
    Us, Zs, T1s, res = [], [], [], {}
    for i, name in enumerate(names):
        P = np.load(Path(args.params_dir) / f"bench_params_{name}.npz")
        U = np.stack([polar(u) for u in P[f"U_{args.tag}"].astype(np.complex128)])
        Z = P[f"Z_{args.tag}"].astype(np.float64)
        rng = np.random.default_rng([args.seed, i])
        gU, gZ = [U], [Z]
        for _ in range(args.group):
            a, b = perturb(U, Z, rng)
            gU.append(a)
            gZ.append(b)
        Us.append(np.stack(gU))
        Zs.append(np.stack(gZ))
        T1s.append(np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t1"].astype(np.float64))
        res[name] = [None] * (args.group + 1)
    np.savez(args.out + ".npz", names=np.array(names), U=np.stack(Us), Z=np.stack(Zs), t1=np.stack(T1s),
             tag=args.tag, seed=args.seed)
    json.dump({"args": vars(args), "per_molecule": res}, open(args.out + ".json", "w"))


if __name__ == "__main__":
    main()
