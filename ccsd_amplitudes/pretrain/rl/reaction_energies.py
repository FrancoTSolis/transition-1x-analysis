#!/usr/bin/env python3
"""Reaction energies and barriers from existing exact LUCJ energies (norb 17-19 R/P/TS triplets), zero compute.

For every reaction with all three geometries (R, P, TS) and every candidate parameter set, computes
    dE_rxn = E(P) - E(R),   dE_barrier = E(TS) - E(R)
and its error against CCSD(T) (pretrain/opt_true/results/ccsd_t_refs.json) and CCSD, in kcal/mol.  Also reports the
spread of the correlation-energy fraction across R/P/TS (non-parallelity).

  python3 -m pretrain.rl.reaction_energies [--out pretrain/opt_true/results/reaction_energies_n17_19.json]
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain" / "opt_true" / "results"
HA2KCAL = 627.5094740631
FILES = ["energy_largeval_base.json", "energy_largeval_pre4_rl4L.json", "energy_largeval_rl4n29.json",
         "energy_largeval_rl4n29_C2H6O.json", "energy_n19_all.json", "energy_n19_rl4n29.json"]
ALIAS = {"rl1": "rl1s", "rl4": "rl4s"}          # the norb 17-18 files used the short names
LAB = {"label": "optimize=True label", "pre1": "pretrained, one shot", "pre4": "pretrained, 4 recycles",
       "rl1s": "small-set RL, one shot", "rl4s": "small-set RL, 4 recycles", "rl4L": "norb 15-18 RL (step 25)",
       "rl4Lf": "norb 15-18 RL (step 50)", "rl4n29": "+ n29 TN RL (step 15)", "rl4n29f": "+ n29 TN RL (step 30)",
       "HF": "Hartree-Fock", "CCSD": "CCSD"}


def load():
    E = defaultdict(dict)
    for f in FILES:
        p = RES / f
        if not p.exists():
            continue
        for n, d in json.load(open(p))["per_molecule"].items():
            for c, v in d.items():
                if not v.get("error") and v.get("E") is not None and np.isfinite(v["E"]):
                    E[n][ALIAS.get(c, c)] = float(v["E"])
    return E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RES / "reaction_energies_n17_19.json"))
    args = ap.parse_args()
    E = load()
    refs = json.load(open(RES / "ccsd_t_refs.json"))
    for n in E:
        if n in refs and "e_ccsd_t" in refs[n]:
            E[n]["HF"], E[n]["CCSD"], E[n]["CCSD(T)"] = refs[n]["e_hf"], refs[n]["e_ccsd"], refs[n]["e_ccsd_t"]
    rx = defaultdict(dict)
    for n in E:
        f, r, g = n.split("_")
        rx[f"{f}_{r}"][g] = n
    full = {k: v for k, v in rx.items() if all(g in v for g in ("R", "P", "TS"))}
    cands = [c for c in LAB if all(c in E[v[g]] for v in full.values() for g in ("R", "P", "TS"))]
    out = {"reactions": {}, "summary": {}}
    for k, v in sorted(full.items()):
        ref = {q: E[v[q]]["CCSD(T)"] for q in ("R", "P", "TS")}
        row = {"ref_ccsd_t": {"rxn": (ref["P"] - ref["R"]) * HA2KCAL, "barrier": (ref["TS"] - ref["R"]) * HA2KCAL}}
        for c in cands:
            e = {q: E[v[q]][c] for q in ("R", "P", "TS")}
            row[c] = {"rxn_err": ((e["P"] - e["R"]) - (ref["P"] - ref["R"])) * HA2KCAL,
                      "barrier_err": ((e["TS"] - e["R"]) - (ref["TS"] - ref["R"])) * HA2KCAL}
        out["reactions"][k] = row
    print(f"{len(full)} complete R/P/TS reactions (norb 17-19); reference CCSD(T); errors in kcal/mol\n")
    print("| candidate | reaction energy MAE | max | barrier MAE | max | mean signed barrier error |")
    print("|:--|--:|--:|--:|--:|--:|")
    for c in cands:
        re_ = np.array([out["reactions"][k][c]["rxn_err"] for k in out["reactions"]])
        be = np.array([out["reactions"][k][c]["barrier_err"] for k in out["reactions"]])
        out["summary"][c] = {"rxn_mae": float(np.abs(re_).mean()), "rxn_max": float(np.abs(re_).max()),
                             "barrier_mae": float(np.abs(be).mean()), "barrier_max": float(np.abs(be).max()),
                             "barrier_mse": float(be.mean())}
        s = out["summary"][c]
        print(f"| {LAB.get(c, c)} | {s['rxn_mae']:.1f} | {s['rxn_max']:.1f} | {s['barrier_mae']:.1f} | "
              f"{s['barrier_max']:.1f} | {s['barrier_mse']:+.1f} |")
    print("\nper reaction, CCSD(T) reaction energy / barrier (kcal/mol):")
    for k, row in out["reactions"].items():
        print(f"  {k}: {row['ref_ccsd_t']['rxn']:+.1f} / {row['ref_ccsd_t']['barrier']:+.1f}")
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"\n-> {args.out}")


if __name__ == "__main__":
    main()
