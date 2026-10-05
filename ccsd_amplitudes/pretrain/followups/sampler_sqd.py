#!/usr/bin/env python3
"""SQD/QSCI stage of the sampler validation: Lin et al.'s noiseless QSCI on the saved sample sets.

Settings (tn_sampler.LIN_SQD): 10^5 samples, postselection only (max_iterations 1, no configuration recovery),
10 batches of 4,000 distinct configurations drawn without replacement by empirical probability
(qiskit_addon_sqd.subsample), alpha and beta strings merged (symmetrize_spin), truncated to max_dim = 4,000 by
marginal count, PySCF selected CI in the product space (qiskit_addon_sqd.fermion.solve_sci; spin_sq = 0.0 as in
her state-vector task lucj_compressed_t2_task_sci.py, --spin-sq none for her quimb task).  The subspaces come from
tn_sampler.sqd_batches (the library's own first-iteration code, seed 0); her returned energy is the minimum over the
batches.  Identical subspaces are diagonalized once.

Sample sets (from pretrain/followups/sampler_validate.py): "mo" (exact, MO basis -> MO integrals), "exact1",
"exact2" (exact, sampled chi basis), "mps<chi>" (MPS, chi basis) -> chi-basis integrals sampled_basis_integrals(Fr).

  OMP_NUM_THREADS=4 ... python -m pretrain.followups.sampler_sqd --names C2H3N_rxn2858_P --sets mo exact1 mps256
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_sampler as S  # noqa: E402

RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "sampler"
LOGS = ROOT / "rl_runs" / "followups" / "sampler"


def subspace_key(ci):
    a, b = ci
    return hashlib.sha1(np.asarray(a, dtype=np.int64).tobytes() + b"|" + np.asarray(b, dtype=np.int64).tobytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--cand", default="rl4n29f")
    ap.add_argument("--sets", nargs="+", default=["mo", "exact1", "exact2", "mps64", "mps128", "mps256", "mps512"])
    ap.add_argument("--n-batches", type=int, default=None, help="diagonalize only the first n batches")
    ap.add_argument("--spin-sq", default="0.0")
    ap.add_argument("--max-dim", type=int, default=S.LIN_SQD["max_dim"])
    ap.add_argument("--samples-per-batch", type=int, default=S.LIN_SQD["samples_per_batch"])
    ap.add_argument("--tag", default="")
    ap.add_argument("--npz-suffix", default="", help="read <name>__<cand><suffix>.npz (e.g. __seeds)")
    ap.add_argument("--all-sets", action="store_true", help="every ints_* key of the npz (overrides --sets)")
    a = ap.parse_args()
    spin_sq = None if a.spin_sq.lower() == "none" else float(a.spin_sq)
    for name in a.names:
        k_tag = f"{name}__{a.cand}"
        z = np.load(LOGS / "samples" / f"{k_tag}{a.npz_suffix}.npz")
        sets = sorted(k[5:] for k in z.files if k.startswith("ints_")) if a.all_sets else a.sets
        d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
        h, g, const = d["one_body"], d["two_body"], float(d["constant"])
        norb, nelec = int(d["norb"]), (int(d["nelec_a"]), int(d["nelec_b"]))
        hF, gF = S.sampled_basis_integrals(h, g, z["Fr"])
        out = RES / f"sqd_{k_tag}{a.tag}.json"
        rec = json.load(open(out)) if out.exists() else {"name": name, "cand": a.cand, "norb": norb,
                                                         "nelec": list(nelec), "sets": {}}
        rec["settings"] = {**S.LIN_SQD, "max_dim": a.max_dim, "samples_per_batch": a.samples_per_batch,
                           "spin_sq": spin_sq, "seed": 0, "solver": "qiskit_addon_sqd.fermion.solve_sci (pyscf SCI)"}
        cache = {}
        for st in sets:
            key = f"ints_{st}"
            if key not in z:
                print(f"[{name}] {st}: no samples ({key})", flush=True)
                continue
            ints = z[key]
            basis = "mo" if st.startswith("mo") else "chi"
            h1, g1 = (h, g) if basis == "mo" else (hF, gF)
            batches = S.sqd_batches(ints, norb, nelec, seed=0, max_dim=a.max_dim,
                                    samples_per_batch=a.samples_per_batch)
            nb = len(batches) if a.n_batches is None else min(a.n_batches, len(batches))
            prev = rec["sets"].get(st, {})
            done = prev.get("batches", [])
            res = {"basis": basis, "unique": S.unique_counts(ints, norb),
                   "dims": [len(x) * len(y) for x, y in batches], "batches": done}
            for b in range(len(done), nb):
                ck = (basis, subspace_key(batches[b]))
                if ck in cache:
                    E, info = cache[ck]
                    info = {**info, "reused": True}
                else:
                    E, info = S.solve_batch(h1, g1, const, batches[b], norb, nelec, spin_sq=spin_sq)
                    cache[ck] = (E, info)
                res["batches"].append({"E": E, **info})
                Es = [x["E"] for x in res["batches"]]
                res.update(E_min=min(Es), E_mean=float(np.mean(Es)), E_std=float(np.std(Es)), n_done=len(Es))
                rec["sets"][st] = res
                RES.mkdir(parents=True, exist_ok=True)
                json.dump(rec, open(out, "w"), indent=1)
                print(f"[{name}] {st:8s} batch {b}: E {E:.8f}  dim {info['dim_a']}x{info['dim_b']}={info['dim']}  "
                      f"S^2 {info['spin_sq']:.2e}  t {info['t']:.1f}s{' (reused)' if info.get('reused') else ''}",
                      flush=True)
            rec["sets"][st] = res
            json.dump(rec, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
