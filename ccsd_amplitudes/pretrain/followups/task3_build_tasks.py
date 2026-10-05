#!/usr/bin/env python3
"""Task 3: energy-task pickles (policy_dump format) for the converged optimize=True labels.

Candidates per molecule (tags):
    labconv        final parameters of the long fit (L-BFGS <= --maxiter, rhf_targets_compressed_conv)
    lab<k>         snapshot at iteration k of the same fit (e.g. lab1000), only if the fit ran past k
    label_ctl      the existing 500-iteration label, re-evaluated as a control (only for --control-names)
The original 500-iteration label is the snapshot at 500 of the same trajectory (bit-identical, checked by
task3_fit_labels.py); its energies are taken from the existing result JSONs.

    pretrain/.train_venv/bin/python3 pretrain/followups/task3_build_tasks.py --names-file pretrain/rl/n29_test24.txt \
        --snap 1000 --control-names C5H5N_rxn3303_P --control-dir rhf_targets_compressed_n29test --out runs_ot/...pkl
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names-file", nargs="+", required=True)
    ap.add_argument("--labels-dir", default="rhf_targets_compressed_conv")
    ap.add_argument("--config", default="square_reg0.005")
    ap.add_argument("--no-final", action="store_true", help="do not add the final (labconv) candidate")
    ap.add_argument("--snap", type=int, nargs="*", default=[], help="snapshot iterations to add as lab<k>")
    ap.add_argument("--min-nit", type=int, default=0,
                    help="add labconv only if the fit ran more than this many iterations (0: always)")
    ap.add_argument("--only-names", nargs="*", default=None)
    ap.add_argument("--control-names", nargs="*", default=[])
    ap.add_argument("--control-dir", default=None, help="dir of the existing labels for --control-names")
    ap.add_argument("--prefix", default="lab",
                    help="tag prefix: final -> <prefix>conv, snapshot k -> <prefix><k> (default labconv / lab500)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    names = []
    for nf in args.names_file:
        names += [ln.strip() for ln in open(ROOT / nf) if ln.strip()]
    names = list(dict.fromkeys(names))
    if args.only_names is not None:
        names = [n for n in names if n in set(args.only_names)]
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    tasks, meta = [], {}
    for n in names:
        lab = np.load(ROOT / args.labels_dir / args.config / f"{n}.npz")
        t1 = np.load(ROOT / "rhf_dataset" / f"{n}.npz")["t1"].astype(np.float64)
        m = meta.setdefault(n, {"resid": {}, "norb": idx[n][0], "nit": int(lab["nit"]),
                                "success": bool(lab["success"])})
        if not args.no_final and int(lab["nit"]) > args.min_nit:
            tasks.append(((n, f"{args.prefix}conv"), n, lab["U_re"] + 1j * lab["U_im"], lab["Z"], t1))
            m["resid"][f"{args.prefix}conv"] = float(lab["resid"])
        its = list(lab["snap_iters"])
        for k in args.snap:
            if k in its:
                j = its.index(k)
                tasks.append(((n, f"{args.prefix}{k}"), n, lab["snap_U_re"][j] + 1j * lab["snap_U_im"][j],
                              lab["snap_Z"][j], t1))
                m["resid"][f"{args.prefix}{k}"] = float(lab["snap_resid"][j])
        if n in args.control_names:
            ref = np.load(ROOT / args.control_dir / args.config / f"{n}.npz")
            tasks.append(((n, "label_ctl"), n, ref["U_re"] + 1j * ref["U_im"], ref["Z"], t1))
            m["resid"]["label_ctl"] = float(ref["resid"])
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    pickle.dump({"tasks": tasks, "meta": meta}, open(out, "wb"))
    tags = sorted({t[0][1] for t in tasks})
    print(f"{len(tasks)} tasks ({', '.join(f'{g}: {sum(t[0][1] == g for t in tasks)}' for g in tags)}) -> {args.out}")


if __name__ == "__main__":
    main()
