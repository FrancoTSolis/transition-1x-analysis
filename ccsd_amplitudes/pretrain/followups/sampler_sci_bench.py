#!/usr/bin/env python3
"""Cost of one PySCF selected-CI sigma build (H c in the fixed space strs x strs) for the norb-29 TN-sampled QSCI
subspaces, truncated to the max_dim most frequent strings (Lin et al.'s ordering), to extrapolate the full batch.

  OMP_NUM_THREADS=2 python -m pretrain.followups.sampler_sci_bench --name C3H5N3_rxn4472_P --set mps256 --dims 100 200 400
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_sampler as S  # noqa: E402

RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "sampler"
LOGS = ROOT / "rl_runs" / "followups" / "sampler"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="C3H5N3_rxn4472_P")
    ap.add_argument("--cand", default="rl4n29f")
    ap.add_argument("--set", default="mps256")
    ap.add_argument("--dims", type=int, nargs="+", default=[100, 200, 400])
    ap.add_argument("--reps", type=int, default=2)
    a = ap.parse_args()
    from pyscf.fci import selected_ci
    z = np.load(LOGS / "samples" / f"{a.name}__{a.cand}.npz")
    d = np.load(ROOT / "rhf_hamiltonians" / f"{a.name}.npz")
    norb, nelec = int(d["norb"]), (int(d["nelec_a"]), int(d["nelec_b"]))
    hF, gF = S.sampled_basis_integrals(d["one_body"], d["two_body"], z["Fr"])
    out = {"name": a.name, "set": a.set, "threads": os.environ["OMP_NUM_THREADS"], "rows": []}
    full = S.sqd_batches(z[f"ints_{a.set}"], norb, nelec, seed=0)[0]
    out["full_strings"] = [len(full[0]), len(full[1])]
    for md in a.dims:
        sa, sb = S.sqd_batches(z[f"ints_{a.set}"], norb, nelec, seed=0, max_dim=md)[0]
        myci = selected_ci.SelectedCI()
        t0 = time.time()
        link = selected_ci._all_linkstr_index((sa, sb), norb, nelec)
        t_link = time.time() - t0
        h2e = myci.absorb_h1e(hF, gF, norb, nelec, 0.5)
        rng = np.random.default_rng(0)
        c = selected_ci._as_SCIvector(rng.normal(size=(len(sa), len(sb))), (sa, sb))
        ts = []
        for _ in range(a.reps):
            t0 = time.time()
            myci.contract_2e(h2e, c, norb, nelec, link)
            ts.append(time.time() - t0)
        t0 = time.time()
        selected_ci.contract_ss(c, norb, nelec)
        t_ss = time.time() - t0
        row = {"max_dim": md, "n_a": len(sa), "n_b": len(sb), "dim": len(sa) * len(sb),
               "dd_inter_a": int(link[1].shape[0]), "t_link": t_link, "t_sigma": min(ts), "t_ss": t_ss}
        out["rows"].append(row)
        print(row, flush=True)
        json.dump(out, open(RES / f"sci_bench_{a.name}__{a.cand}_{a.set}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
