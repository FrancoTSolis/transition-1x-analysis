#!/usr/bin/env python3
"""Repeated sample sets for the QSCI statistics of the sampler validation (norb 15-16).

QSCI energies from 10^5 samples fluctuate from sample set to sample set (which rare alpha/beta strings happen to be
drawn).  To compare MPS with exact sampling beyond one draw, this script draws R independent sets of 10^5 samples
from the exact chi-basis state and from the MPS at each chi (fresh seeds, disjoint from sampler_validate.py), and
saves them to rl_runs/followups/sampler/samples/<name>__<cand>__seeds.npz (keys ints_exact_s<r>, ints_mps<chi>_s<r>)
for pretrain/followups/sampler_sqd.py --npz-suffix __seeds.

  CUDA_VISIBLE_DEVICES=5 OMP_NUM_THREADS=2 ... python -m pretrain.followups.sampler_seeds --names ... --chis 64 256 -R 6
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_sampler as S  # noqa: E402
from pretrain.followups.sampler_validate import basis_cache_for, load_params  # noqa: E402

LOGS = ROOT / "rl_runs" / "followups" / "sampler"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--cand", default="rl4n29f")
    ap.add_argument("--tasks", default="runs_ot/energy_tasks/smallval_rl4n29.pkl")
    ap.add_argument("--chis", type=int, nargs="*", default=[64, 256])
    ap.add_argument("--mo", action="store_true", help="also draw exact MO-basis sets (ints_mo_s<r>, MO integrals)")
    ap.add_argument("-R", type=int, default=6)
    ap.add_argument("--r0", type=int, default=0, help="set indices r0+1 .. r0+R (fresh seeds for extra sets)")
    ap.add_argument("--suffix", default="__seeds", help="output npz <name>__<cand><suffix>.npz")
    ap.add_argument("--no-exact", action="store_true")
    ap.add_argument("--shots", type=int, default=S.LIN_SHOTS)
    ap.add_argument("--max-mem-gb", type=float, default=4.0)
    a = ap.parse_args()
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
    for name in a.names:
        t0 = time.time()
        d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
        h, g, const = d["one_body"], d["two_body"], float(d["constant"])
        norb, nelec = int(d["norb"]), (int(d["nelec_a"]), int(d["nelec_b"]))
        U, Z, t1 = load_params(ROOT / a.tasks, name, a.cand)
        from pretrain.rl.tn_energy import LUCJEnergyTN
        ev = LUCJEnergyTN(h, g, const, norb, nelec, max_bond=max(a.chis or [64]), device="cuda", name=name,
                          basis_cache=basis_cache_for(name), block2_threads=1,
                          scratch=tempfile.mkdtemp(prefix="smps_", dir=str(LOGS)), stack_mem=1 << 30)
        Fr = S.sampled_basis(ev, t1)
        hF, gF = S.sampled_basis_integrals(h, g, Fr)
        save = {"Fr": Fr, "occ_mask": np.asarray(ev.occ), "norb": norb, "nelec": np.asarray(nelec)}
        strs = S.ci_strings(norb, nelec[0])
        if a.mo:
            eng = S.LUCJStateGPU(h, g, const, norb, nelec, max_mem_gb=a.max_mem_gb)
            P = eng.state(U, Z, t1)
            eng.release()
            for r in range(a.r0 + 1, a.r0 + a.R + 1):
                ia, ib = S.sample_ci_matrix(P, a.shots, seed=70_000 + r)
                save[f"ints_mo_s{r}"] = strs[ia] | (strs[ib] << norb)
            del P
            torch.cuda.empty_cache()
        eng = S.LUCJStateGPU(hF, gF, const, norb, nelec, basis=Fr, max_mem_gb=a.max_mem_gb)
        P = eng.state(U, Z, t1)
        eng.release()
        for r in range(a.r0 + 1, a.r0 + a.R + 1):
            if a.no_exact:
                break
            ia, ib = S.sample_ci_matrix(P, a.shots, seed=50_000 + r)
            save[f"ints_exact_s{r}"] = strs[ia] | (strs[ib] << norb)
        del P
        torch.cuda.empty_cache()
        for chi in a.chis:
            mps, _, info = S.build_mps(ev, U, Z, t1, chi)
            smp = S.MPSSampler(mps)
            for r in range(a.r0 + 1, a.r0 + a.R + 1):
                occ, _ = smp.sample(a.shots, seed=60_000 + 100 * chi + r, return_logp=False)
                save[f"ints_mps{chi}_s{r}"] = S.occ_to_ints(occ)
            del mps, smp
            torch.cuda.empty_cache()
        (LOGS / "samples").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(LOGS / "samples" / f"{name}__{a.cand}{a.suffix}.npz", **save)
        print(f"[{name}] {a.R} sets each of exact + MPS chi {a.chis}: {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
