#!/usr/bin/env python3
"""End-to-end TN-sampled QSCI demonstration at norb 29 (58 qubits): dm-engine MPS -> 10^5 direct samples ->
Lin et al.'s SQD subsampling (10 x 4,000) -> subspace dimensions [-> QSCI diagonalization: pretrain/followups/
sampler_sqd.py on the saved samples].  No exact reference exists at this size, so the sample distribution at chi is
checked against a larger-chi MPS of the same state (Monte-Carlo estimators from the larger-chi samples):
    TVD(p_chi, p_ref) = 1/2 E_{x~p_ref} |1 - p_chi(x)/p_ref(x)|,   |<psi_ref|psi_chi>|^2 = |E_{x~p_ref}[psi_chi/psi_ref]|^2.

  CUDA_VISIBLE_DEVICES=5 OMP_NUM_THREADS=4 ... python -m pretrain.followups.sampler_n29 --name C3H5N3_rxn4472_P \\
      --chis 128 256 512 --ref-chi 512
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import sys
import tempfile
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl import tn_sampler as S  # noqa: E402
from pretrain.followups.sampler_validate import basis_cache_for, load_params, mc_mean  # noqa: E402

RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "sampler"
LOGS = ROOT / "rl_runs" / "followups" / "sampler"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="C3H5N3_rxn4472_P")
    ap.add_argument("--cand", default="rl4n29f")
    ap.add_argument("--tasks", default="runs_ot/energy_tasks/n29test_rl4n29.pkl")
    ap.add_argument("--chis", type=int, nargs="+", default=[128, 256, 512])
    ap.add_argument("--ref-chi", type=int, default=512)
    ap.add_argument("--shots", type=int, default=S.LIN_SHOTS)
    ap.add_argument("--block2-threads", type=int, default=4)
    ap.add_argument("--no-energy", action="store_true")
    ap.add_argument("--throughput-shots", type=int, default=1_000_000, help="extra timing run (not saved), 0: off")
    a = ap.parse_args()
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "4")))
    name = a.name
    k_tag = f"{name}__{a.cand}"
    out_json = RES / f"n29_{k_tag}.json"
    out_npz = LOGS / "samples" / f"{k_tag}.npz"
    d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
    h, g, const = d["one_body"], d["two_body"], float(d["constant"])
    norb, nelec = int(d["norb"]), (int(d["nelec_a"]), int(d["nelec_b"]))
    U, Z, t1 = load_params(ROOT / a.tasks, name, a.cand)
    refs = json.load(open(ROOT / "pretrain" / "opt_true" / "results" / "ccsd_t_refs.json")).get(name, {})
    prev = json.load(open(ROOT / "pretrain" / "opt_true" / "results" / "energy_n29test_rl4n29_chi256.json"))
    e_v1 = prev["per_molecule"].get(name, {}).get(a.cand, {}).get("E")
    rec = {"name": name, "cand": a.cand, "norb": norb, "nelec": list(nelec), "e_hf": float(d["e_hf"]),
           "e_ccsd": float(d["e_ccsd"]), "e_ccsd_t": refs.get("e_ccsd_t"), "e_tn_v1_chi256_m1.5": e_v1,
           "shots": a.shots, "host": socket.gethostname(), "gpu": torch.cuda.get_device_name(), "lin_sqd": S.LIN_SQD,
           "mps": {}}
    from pretrain.rl.tn_energy import LUCJEnergyTN
    ev = LUCJEnergyTN(h, g, const, norb, nelec, max_bond=max(a.chis), device="cuda", name=name,
                      basis_cache=basis_cache_for(name), block2_threads=a.block2_threads,
                      scratch=tempfile.mkdtemp(prefix="smp29_", dir=str(LOGS)), stack_mem=2 << 30)
    Fr = S.sampled_basis(ev, t1)
    rec["tn_settings"] = ev.settings()
    save = {"Fr": Fr, "occ_mask": np.asarray(ev.occ), "norb": norb, "nelec": np.asarray(nelec)}
    hf_int = int(S.occ_to_ints(np.array([[3 if o else 0 for o in ev.occ]], dtype=np.uint8))[0])
    samplers, occs = {}, {}
    for chi in sorted(a.chis, key=lambda c: c != a.ref_chi):         # reference chi first
        t0 = time.time()
        mps, Fr2, info = S.build_mps(ev, U, Z, t1, chi)
        assert np.allclose(Fr2, Fr)
        r = {kk: info[kk] for kk in ("discarded_sum", "discarded_max", "n_trunc", "max_bond", "bond_dims", "t_state")}
        if not a.no_energy:
            E, einfo = S.mps_energy(ev, mps, Fr)
            r.update(E=E, corr_frac=(E - rec["e_hf"]) / (rec["e_ccsd"] - rec["e_hf"]), t_energy=einfo["t_mpo"] + einfo["t_expect"])
        t1_ = time.time()
        smp = S.MPSSampler(mps)
        torch.cuda.synchronize()
        t2_ = time.time()
        occ, lp = smp.sample(a.shots, seed=2000 + chi)
        torch.cuda.synchronize()
        t3_ = time.time()
        ints = S.occ_to_ints(occ)
        save[f"ints_mps{chi}"] = ints
        vals, cnt = np.unique(ints, return_counts=True)
        top = np.argsort(-cnt)[:10]
        batches = S.sqd_batches(ints, norb, nelec, seed=0)
        r.update(t_canonicalize=t2_ - t1_, t_sample=t3_ - t2_, isometry_err=smp.isometry_err,
                 unique=S.unique_counts(ints, norb), freq_top10=(cnt[top] / a.shots).tolist(),
                 top1_is_hf_s=bool(vals[top[0]] == hf_int), p_hf_s=float(np.exp(lp[ints == hf_int][0])) if
                 (ints == hf_int).any() else None,
                 sqd_dims=[len(x) * len(y) for x, y in batches],
                 sqd_dims_ab=[[len(x), len(y)] for x, y in batches])
        if a.throughput_shots:
            torch.cuda.synchronize()
            tt = time.time()
            smp.sample(a.throughput_shots, seed=3000 + chi, return_logp=False)
            torch.cuda.synchronize()
            r.update(throughput_shots=a.throughput_shots, t_sample_throughput=time.time() - tt)
        samplers[chi], occs[chi] = smp, occ
        if chi != a.ref_chi and a.ref_chi in samplers:
            ref, oref = samplers[a.ref_chi], occs[a.ref_chi]
            sg = S.ffsim_sign(oref)                                   # same basis/convention: sign cancels; kept
            p_ref = ref.amplitudes(oref) * sg
            p_chi = smp.amplitudes(oref) * sg
            ov, ov_se = mc_mean(p_chi / p_ref)
            tv, tv_se = mc_mean(0.5 * np.abs(1 - np.abs(p_chi) ** 2 / np.abs(p_ref) ** 2))
            r.update(vs_ref_chi=a.ref_chi, fidelity_vs_ref=abs(ov) ** 2, overlap_se=ov_se, tvd_vs_ref=tv, tvd_vs_ref_se=tv_se)
        rec["mps"][str(chi)] = r
        print(f"[{name}] chi {chi}: E {r.get('E', float('nan')):.6f} corr {100 * r.get('corr_frac', float('nan')):.2f}% "
              f"disc {r['discarded_sum']:.2e} t_state {r['t_state']:.0f}s sample {r['t_sample']:.2f}s iso "
              f"{smp.isometry_err:.1e} unique {r['unique']}  sqd dims {r['sqd_dims']}  "
              + (f"F(vs {a.ref_chi}) {r['fidelity_vs_ref']:.5f} TVD {r['tvd_vs_ref']:.4f}" if "fidelity_vs_ref" in r else ""),
              flush=True)
        del mps
        torch.cuda.empty_cache()
        RES.mkdir(parents=True, exist_ok=True)
        json.dump(rec, open(out_json, "w"), indent=1)
        (LOGS / "samples").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_npz, **save)
    print(f"-> {out_json}", flush=True)


if __name__ == "__main__":
    main()
