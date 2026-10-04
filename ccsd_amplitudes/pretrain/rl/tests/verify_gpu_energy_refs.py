#!/usr/bin/env python3
"""[verifier] Stored-reference validation of pretrain/rl/gpu_energy.py on an independent random task selection.

norb 15-16: a seeded random subset of the 390 tasks (all candidate types mixed), complex64 with engine reuse
(set_npz) and complex128 with fresh engines; reports |dE| vs the stored CPU energies AND vs the engine's own
complex128 value (c8 precision without the stored-reference noise).  Also checks energies() == energy() and
set_npz-vs-fresh-engine bit identity.

  CUDA_VISIBLE_DEVICES=3 python -m pretrain.rl.tests.verify_gpu_energy_refs --n 48 --seed 2026 --out X.json
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402

from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac  # noqa: E402

TAGS = ["baselines", "all1", "all4", "all8", "swap"]


def load_small():
    out = []
    for tag in TAGS:
        tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{tag}.pkl", "rb"))["tasks"]
        ref = json.load(open(ROOT / "pretrain/opt_true/results" / f"energy_small_{tag}.json"))["per_molecule"]
        for (name, cand), _, U, Z, t1 in tasks:
            r = ref.get(name, {}).get(cand)
            if r is None:
                continue
            out.append(dict(tag=tag, name=name, cand=cand, U=U, Z=Z, t1=t1, E_ref=float(r["E"]),
                            cf_ref=float(r["corr_frac"])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=48)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--include", nargs="*", default=[], help="tag:name:cand tasks always included")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rows = load_small()
    rng = np.random.default_rng(args.seed)
    idx = set(rng.choice(len(rows), size=args.n, replace=False).tolist())
    for inc in args.include:
        tag, name, cand = inc.split(":")
        idx |= {i for i, r in enumerate(rows) if (r["tag"], r["name"], r["cand"]) == (tag, name, cand)}
    sel = [rows[i] for i in sorted(idx)]
    sel.sort(key=lambda r: r["name"])
    print(f"{len(sel)} tasks, {len(set(r['name'] for r in sel))} molecules, cands:",
          sorted(set(r["cand"] for r in sel)), flush=True)

    # ---- complex128, fresh engine per molecule
    t_all = time.time()
    by = {}
    for r in sel:
        by.setdefault(r["name"], []).append(r)
    for name, rs in by.items():
        eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex128)
        for r in rs:
            E = eng.energy(r["U"], r["Z"], t1=r["t1"])
            r["E16"] = E
            r["norb"], r["k"] = eng.norb, eng.k
            r["e_hf"], r["e_ccsd"] = eng.e_hf, eng.e_ccsd
        eng.release()
        del eng
        torch.cuda.empty_cache()
    print(f"c16 done ({time.time() - t_all:.0f}s)", flush=True)

    # ---- complex64, one engine per (norb, k) reused through set_npz
    engs = {}
    for name, rs in by.items():
        key = (rs[0]["norb"], rs[0]["k"])
        if key in engs:
            eng = engs[key].set_npz(name)
        else:
            for e in engs.values():
                e.release()
            engs.clear()
            eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex64)
            engs[key] = eng
        for r in rs:
            E = eng.energy(r["U"], r["Z"], t1=r["t1"])
            r["E8"] = E
            r["t8"] = eng.timing["total"]
    # bit-identity of set_npz vs a fresh engine and energies() vs energy() on the last molecule
    last = list(by)[-1]
    rs = by[last]
    eng = engs[(rs[0]["norb"], rs[0]["k"])]
    Es_batch = eng.energies([(r["U"], r["Z"], r["t1"]) for r in rs])
    for e in engs.values():
        e.release()
    engs.clear()
    torch.cuda.empty_cache()
    fresh = LUCJEnergyGPU.from_npz(last, dtype=torch.complex64)
    Es_fresh = [fresh.energy(r["U"], r["Z"], t1=r["t1"]) for r in rs]
    fresh.release()
    ident = dict(name=last, reuse=[r["E8"] for r in rs], batch=Es_batch, fresh=Es_fresh,
                 batch_eq_single=all(a == b for a, b in zip(Es_batch, [r["E8"] for r in rs])),
                 fresh_eq_reuse=all(a == b for a, b in zip(Es_fresh, [r["E8"] for r in rs])))
    print("identity checks:", json.dumps({k: v for k, v in ident.items() if "eq" in k or k == "name"}), flush=True)

    out = []
    for r in sel:
        d = dict(tag=r["tag"], name=r["name"], cand=r["cand"], norb=r["norb"], k=r["k"], E_ref=r["E_ref"],
                 E16=r["E16"], E8=r["E8"], dE16=r["E16"] - r["E_ref"], dE8=r["E8"] - r["E_ref"],
                 dE8_vs16=r["E8"] - r["E16"],
                 dcf8=corr_frac(r["E8"], r["e_hf"], r["e_ccsd"]) - r["cf_ref"],
                 dcf16=corr_frac(r["E16"], r["e_hf"], r["e_ccsd"]) - r["cf_ref"],
                 dcf8_vs16=corr_frac(r["E8"], r["e_hf"], r["e_ccsd"]) - corr_frac(r["E16"], r["e_hf"], r["e_ccsd"]),
                 cf16=corr_frac(r["E16"], r["e_hf"], r["e_ccsd"]), t8=r["t8"])
        out.append(d)
        print(f"  {d['name']:18s} {d['cand']:22s} ({d['norb']},{d['k']})  corr% {100*d['cf16']:7.3f}  "
              f"dE16 {d['dE16']:+.2e}  dE8 {d['dE8']:+.2e}  dE8-16 {d['dE8_vs16']:+.2e}  dcf8 {d['dcf8']:+.1e}"
              f"  t8 {d['t8']:.2f}s", flush=True)
    summ = dict(n=len(out),
                max_abs_dE8=max(abs(d["dE8"]) for d in out), max_abs_dcf8=max(abs(d["dcf8"]) for d in out),
                max_abs_dE16=max(abs(d["dE16"]) for d in out), max_abs_dcf16=max(abs(d["dcf16"]) for d in out),
                max_abs_dE8_vs16=max(abs(d["dE8_vs16"]) for d in out),
                max_abs_dcf8_vs16=max(abs(d["dcf8_vs16"]) for d in out),
                mean_dE16=float(np.mean([d["dE16"] for d in out])))
    for shape in sorted(set((d["norb"], d["k"]) for d in out)):
        dd = [d for d in out if (d["norb"], d["k"]) == shape]
        summ[f"shape_{shape[0]}_{shape[1]}"] = dict(n=len(dd), max_abs_dE8=max(abs(d["dE8"]) for d in dd),
                                                    max_abs_dE16=max(abs(d["dE16"]) for d in dd),
                                                    mean_dE16=float(np.mean([d["dE16"] for d in dd])),
                                                    max_abs_dE8_vs16=max(abs(d["dE8_vs16"]) for d in dd))
    print(json.dumps(summ, indent=1))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(summary=summ, identity=ident, rows=out, args=vars(args)), open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
