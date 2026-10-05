#!/usr/bin/env python3
"""CCSD(T) reference energies in the same frozen-core STO-3G active space as rhf_hamiltonians/ (cached JSON).

Rebuilds the RHF of pretrain/rl/hamiltonian.py (same geometry, basis, frozen core), runs CCSD + (T), and checks that
E_HF and E_CCSD reproduce the stored Hamiltonian file.  Results are merged into the cache file:
    pretrain/opt_true/results/ccsd_t_refs.json : {name: {e_hf, e_ccsd, e_ccsd_t, t_s}}

  python3 -m pretrain.rl.ccsd_t_refs --names-file gauge_study/names_norb19.txt [--n-procs 8]
"""
from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "pretrain" / "opt_true" / "results" / "ccsd_t_refs.json"


def one(name):
    os.environ["OMP_NUM_THREADS"] = "1"
    from pyscf import cc, gto, scf
    from pretrain.rl.hamiltonian import n_frozen_core, read_xyz
    t0 = time.time()
    xyz = read_xyz(ROOT / "jobs" / name / f"{name}.xyz")
    mol = gto.M(atom=xyz, basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf, frozen=n_frozen_core(xyz)).run()
    et = mycc.ccsd_t()
    h = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
    d_hf, d_cc = float(mf.e_tot) - float(h["e_hf"]), float(mycc.e_tot) - float(h["e_ccsd"])
    if abs(d_hf) > 1e-6 or abs(d_cc) > 1e-5:
        return name, {"error": f"mismatch with rhf_hamiltonians: dE_HF {d_hf:.2e}, dE_CCSD {d_cc:.2e}"}
    return name, {"e_hf": float(mf.e_tot), "e_ccsd": float(mycc.e_tot), "e_ccsd_t": float(mycc.e_tot + et),
                  "ccsd_converged": bool(mycc.converged), "t_s": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names-file", action="append", required=True)
    ap.add_argument("--n-procs", type=int, default=8)
    args = ap.parse_args()
    names = []
    for f in args.names_file:
        names += [ln.strip() for ln in open(ROOT / f) if ln.strip()]
    cache = json.load(open(CACHE)) if CACHE.exists() else {}
    todo = [n for n in dict.fromkeys(names) if n not in cache or "error" in cache[n]]
    print(f"{len(names)} names, {len(todo)} to compute", flush=True)
    with Pool(args.n_procs) as pool:
        for name, r in pool.imap_unordered(one, todo):
            cache[name] = r
            print(name, r, flush=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    json.dump(cache, open(CACHE, "w"), indent=1, sort_keys=True)
    print(f"-> {CACHE} ({len(cache)} molecules)")


if __name__ == "__main__":
    main()
