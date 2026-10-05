#!/usr/bin/env python3
"""Validation of the MPS sampler (pretrain/rl/tn_sampler.py) against exact samples at norb 15-16.

For each molecule and candidate (U, Z, t1 from an energy-task pickle):
  1. exact LUCJ state on the GPU in the MO basis (pretrain.rl.gpu_energy) -> E_MO, 10^5 MO-basis samples (for the
     MO-vs-chi QSCI comparison of the SQD stage);
  2. exact state in the sampled basis chi = MO @ Fr, Fr = final(t1) @ S (tn_sampler.LUCJStateGPU) -> E_chi (must equal
     E_MO), two independent sets of 10^5 exact chi-basis samples, the exact top-K configurations;
  3. for every chi: dm-engine MPS (LUCJEnergyTN, Boys/Fiedler basis), block2 energy, 10^5 MPS samples; metrics:
       * TVD on the top-K configurations (all configurations outside the top K pooled into one bin), exact (MPS
         probabilities from contracted amplitudes vs exact probabilities: no sampling noise) and from the samples
         (empirical MPS frequencies vs exact probabilities; the floor = an exact sample set of the same size);
       * full TVD and fidelity |<psi_exact|psi_mps>|^2, unbiased Monte-Carlo estimates from exact samples
         (TVD = 1/2 E_{x~p_e}|1 - p_m/p_e|, <psi_m|psi_e> = E_{x~p_e}[conj(psi_m/psi_e)]) and, as a check of the
         sampler, the overlap estimated from the MPS samples (E_{x~p_m}[psi_e/psi_m]);
       * sampler self-consistency: logp of each draw vs log|amplitude|^2;
       * unique full configurations / alpha / beta / alpha-union-beta strings among 10^5 samples, and the subspace
         dimensions of Lin et al.'s 10 x 4,000-sample SQD batches (tn_sampler.sqd_batches; no diagonalization).
  Samples are saved (npz) for the SQD stage (pretrain/followups/sampler_sqd.py).

  CUDA_VISIBLE_DEVICES=5 OMP_NUM_THREADS=4 ... python -m pretrain.followups.sampler_validate \\
      --names C2H3N_rxn2858_P C2H3N_rxn2858_TS C2H4O_rxn0724_P --cand rl4n29f --chis 64 128 256 512
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

RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "sampler"
LOGS = ROOT / "rl_runs" / "followups" / "sampler"
SHARED_BASIS_CACHE = ROOT / "rl_runs" / "tn_basis_cache"
OWN_BASIS_CACHE = LOGS / "basis_cache"
TOPK = (1, 10, 100, 1000, 10000)


def load_params(tasks_path, name, cand):
    import pickle
    for key, nm, U, Z, t1 in pickle.load(open(tasks_path, "rb"))["tasks"]:
        if nm == name and key[1] == cand:
            return np.asarray(U), np.asarray(Z), None if t1 is None else np.asarray(t1)
    raise KeyError((name, cand))


def basis_cache_for(name):
    if (SHARED_BASIS_CACHE / f"{name}_boys_fiedler.npz").exists():
        return str(SHARED_BASIS_CACHE)            # same basis as the TN reward workers (read only)
    OWN_BASIS_CACHE.mkdir(parents=True, exist_ok=True)
    return str(OWN_BASIS_CACHE)


def tvd_top(p_ref, p_other, K):
    """TVD of the coarse-grained distributions: the first K entries (top-K of p_ref) + one pooled 'rest' bin."""
    a, b = p_ref[:K], p_other[:K]
    return float(0.5 * (np.abs(a - b).sum() + abs((1 - a.sum()) - (1 - b.sum()))))


def freq_on(ints_samples, keys):
    """Empirical frequencies of the configurations `keys` (int64) among the samples."""
    order = np.argsort(keys)
    sk = keys[order]
    pos = np.searchsorted(sk, ints_samples)
    pos = np.clip(pos, 0, len(sk) - 1)
    hit = sk[pos] == ints_samples
    cnt = np.bincount(pos[hit], minlength=len(sk))
    out = np.empty(len(keys))
    out[order] = cnt / len(ints_samples)
    return out


def mc_mean(x):
    """Sample mean and its standard error sqrt(mean|x - mean|^2 / N) (complex or real)."""
    x = np.asarray(x)
    m = x.mean()
    se = float(np.sqrt(np.mean(np.abs(x - m) ** 2) / len(x)))
    return (complex(m) if np.iscomplexobj(x) else float(m)), se


def run_molecule(a, name):
    k_tag = f"{name}__{a.cand}"
    out_json = RES / f"validate_{k_tag}.json"
    out_npz = LOGS / "samples" / f"{k_tag}.npz"
    d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
    h, g, const = d["one_body"], d["two_body"], float(d["constant"])
    norb, nelec = int(d["norb"]), (int(d["nelec_a"]), int(d["nelec_b"]))
    k = nelec[0]
    U, Z, t1 = load_params(ROOT / a.tasks, name, a.cand)
    refs = json.load(open(ROOT / "pretrain" / "opt_true" / "results" / "ccsd_t_refs.json")).get(name, {})
    rec = {"name": name, "cand": a.cand, "norb": norb, "nelec": list(nelec), "e_hf": float(d["e_hf"]),
           "e_ccsd": float(d["e_ccsd"]), "e_ccsd_t": refs.get("e_ccsd_t"), "shots": a.shots, "host": socket.gethostname(),
           "gpu": torch.cuda.get_device_name(), "tasks": a.tasks, "lin_sqd": S.LIN_SQD, "chis": a.chis}
    from pretrain.rl.tn_energy import LUCJEnergyTN
    ev = LUCJEnergyTN(h, g, const, norb, nelec, max_bond=max(a.chis), device="cuda", name=name,
                      basis_cache=basis_cache_for(name), block2_threads=a.block2_threads,
                      scratch=tempfile.mkdtemp(prefix="smpv_", dir=str(LOGS)), stack_mem=2 << 30)
    Fr = S.sampled_basis(ev, t1)
    hF, gF = S.sampled_basis_integrals(h, g, Fr)
    rec["tn_settings"] = ev.settings()
    rec["occ_mask"] = [bool(x) for x in ev.occ]
    strs = S.ci_strings(norb, k)
    dim = len(strs)
    save = {"Fr": Fr, "occ_mask": np.asarray(ev.occ), "norb": norb, "nelec": np.asarray(nelec)}

    # ---- 1. exact, MO basis
    t0 = time.time()
    eng = S.LUCJStateGPU(h, g, const, norb, nelec, max_mem_gb=a.max_mem_gb)
    E_mo = eng.energy(U, Z, t1)
    P = eng.state(U, Z, t1)
    ia, ib = S.sample_ci_matrix(P, a.shots, seed=100)
    save["ints_mo"] = strs[ia] | (strs[ib] << norb)
    pf = (P.real.float() ** 2 + P.imag.float() ** 2).reshape(-1)
    nrm = float(pf.double().sum())
    v, _ = torch.topk(pf, 10)
    rec["exact_mo"] = {"E": E_mo, "corr_frac": (E_mo - rec["e_hf"]) / (rec["e_ccsd"] - rec["e_hf"]),
                       "p_hf": float(pf[0]) / nrm, "p_top10": (v.double().cpu().numpy() / nrm).tolist(),
                       "unique": S.unique_counts(save["ints_mo"], norb), "t": time.time() - t0}
    eng.release()
    del eng, P, pf
    torch.cuda.empty_cache()

    # ---- 2. exact, sampled (chi) basis
    t0 = time.time()
    engF = S.LUCJStateGPU(hF, gF, const, norb, nelec, basis=Fr, max_mem_gb=a.max_mem_gb)
    E_chi = engF.energy(U, Z, t1)
    P = engF.state(U, Z, t1)
    engF.release()
    del engF
    torch.cuda.empty_cache()
    pf = (P.real.float() ** 2 + P.imag.float() ** 2).reshape(-1)
    nrm = float(pf.double().sum())
    sets = {}
    for s_ in (1, 2):
        ia, ib = S.sample_ci_matrix(P, a.shots, seed=s_)
        sets[s_] = (ia, ib, strs[ia] | (strs[ib] << norb))
        save[f"ints_exact{s_}"] = sets[s_][2]
    K = max(TOPK)
    pv, pidx = torch.topk(pf, K)
    p_top = pv.double().cpu().numpy() / nrm
    pidx = pidx.cpu().numpy()
    top_a, top_b = strs[pidx // dim], strs[pidx % dim]
    top_ints = top_a | (top_b << norb)
    occ_top = S.strings_to_occ(top_a, top_b, norb)
    hf_int = int(S.occ_to_ints(np.array([[3 if o else 0 for o in ev.occ]], dtype=np.uint8))[0])
    f_e2 = freq_on(sets[2][2], top_ints)
    rec["exact_chi"] = {"E": E_chi, "E_minus_E_mo": E_chi - E_mo, "norm": nrm, "p_top10": p_top[:10].tolist(),
                        "top1_is_hf_s": bool(top_ints[0] == hf_int), "mass_topK": {str(kk): float(p_top[:kk].sum()) for kk in TOPK},
                        "unique": {str(s_): S.unique_counts(sets[s_][2], norb) for s_ in (1, 2)},
                        "tvd_topK_sample_floor": {str(kk): tvd_top(p_top, f_e2, kk) for kk in TOPK},
                        "sqd_dims": {str(s_): [len(x) * len(y) for x, y in S.sqd_batches(sets[s_][2], norb, nelec, seed=0)]
                                     for s_ in (1, 2)},
                        "t": time.time() - t0}
    rec["exact_mo"]["sqd_dims"] = [len(x) * len(y) for x, y in S.sqd_batches(save["ints_mo"], norb, nelec, seed=0)]
    print(f"[{name}] E_MO {E_mo:.8f}  E_chi - E_MO {E_chi - E_mo:+.1e}  p_HF(MO) {rec['exact_mo']['p_hf']:.4f}  "
          f"p_top1(chi) {p_top[0]:.4f} (HF_S: {rec['exact_chi']['top1_is_hf_s']})  unique(chi) "
          f"{rec['exact_chi']['unique']['1']}  unique(MO) {rec['exact_mo']['unique']}", flush=True)
    # exact amplitudes of exact sample set 1 (for the MC estimators)
    ia1, ib1, ints1 = sets[1]
    ae1 = P[torch.as_tensor(ia1, device=P.device), torch.as_tensor(ib1, device=P.device)].cpu().numpy().astype(
        np.complex128) / np.sqrt(nrm)
    occ1 = S.strings_to_occ(strs[ia1], strs[ib1], norb)
    sgn1 = S.ffsim_sign(occ1)

    # ---- 3. MPS at every chi
    rec["mps"] = {}
    for chi in a.chis:
        t0 = time.time()
        mps, Fr2, info = S.build_mps(ev, U, Z, t1, chi)
        assert np.allclose(Fr2, Fr)
        E_tn, einfo = S.mps_energy(ev, mps, Fr)
        t1_ = time.time()
        smp = S.MPSSampler(mps)
        torch.cuda.synchronize()
        t2_ = time.time()
        occ_m, lp_m = smp.sample(a.shots, seed=1000 + chi)
        torch.cuda.synchronize()
        t3_ = time.time()
        ints_m = S.occ_to_ints(occ_m)
        save[f"ints_mps{chi}"] = ints_m
        pm_top = smp.probabilities(occ_top)
        t4_ = time.time()
        am1 = smp.amplitudes(occ1) * sgn1                    # MPS amplitudes (ffsim convention) of exact samples
        t5_ = time.time()
        r1 = np.conj(am1) / np.conj(ae1)
        ov1, ov1_se = mc_mean(r1)
        tvd_full, tvd_full_se = mc_mean(0.5 * np.abs(1.0 - np.abs(am1) ** 2 / np.abs(ae1) ** 2))
        am_m = smp.amplitudes(occ_m) * S.ffsim_sign(occ_m)
        a_m, b_m = S.occ_to_strings(occ_m)
        iam, ibm = S.ci_addresses(norb, k, a_m), S.ci_addresses(norb, k, b_m)
        ae_m = P[torch.as_tensor(iam, device=P.device), torch.as_tensor(ibm, device=P.device)].cpu().numpy().astype(
            np.complex128) / np.sqrt(nrm)
        ov2, ov2_se = mc_mean(ae_m / am_m)
        dlogp = np.abs(lp_m - np.log(np.abs(am_m) ** 2))
        f_m = freq_on(ints_m, top_ints)
        r = {"E": E_tn, "E_minus_exact_mHa": 1e3 * (E_tn - E_mo),
             "corr_frac": (E_tn - rec["e_hf"]) / (rec["e_ccsd"] - rec["e_hf"]),
             **{kk: info[kk] for kk in ("discarded_sum", "discarded_max", "n_trunc", "max_bond", "bond_dims")},
             "t_state": info["t_state"], "t_energy": t1_ - t0 - info["t_state"], "t_canonicalize": t2_ - t1_,
             "t_sample": t3_ - t2_, "t_amp_topK": t4_ - t3_, "t_amp_1e5": t5_ - t4_,
             "isometry_err": smp.isometry_err, "norm2_before": smp.norm2,
             "fidelity_from_exact_samples": abs(ov1) ** 2, "overlap_abs_se": ov1_se,
             "fidelity_from_mps_samples": abs(ov2) ** 2, "overlap2_abs_se": ov2_se,
             "tvd_full_mc": tvd_full, "tvd_full_mc_se": tvd_full_se,
             "tvd_topK_exact": {str(kk): tvd_top(p_top, pm_top, kk) for kk in TOPK},
             "tvd_topK_samples_vs_exact_p": {str(kk): tvd_top(p_top, f_m, kk) for kk in TOPK},
             "tvd_topK_samples_vs_exact_samples": {str(kk): tvd_top(f_e2, f_m, kk) for kk in TOPK},
             "mass_topK_mps": {str(kk): float(pm_top[:kk].sum()) for kk in TOPK},
             "max_rel_err_p_top100": float(np.max(np.abs(pm_top[:100] / p_top[:100] - 1))),
             "logp_selfconsistency_max": float(dlogp.max()), "logp_selfconsistency_median": float(np.median(dlogp)),
             "unique": S.unique_counts(ints_m, norb),
             "sqd_dims": [len(x) * len(y) for x, y in S.sqd_batches(ints_m, norb, nelec, seed=0)]}
        rec["mps"][str(chi)] = r
        print(f"[{name}] chi {chi:4d}: E-E_ex {r['E_minus_exact_mHa']:+.3f} mHa  disc {r['discarded_sum']:.2e}  "
              f"F(ex) {r['fidelity_from_exact_samples']:.6f} F(mps) {r['fidelity_from_mps_samples']:.6f}  "
              f"TVD full {tvd_full:.4f}+-{tvd_full_se:.4f}  TVD top1e3 exact {r['tvd_topK_exact']['1000']:.4f} "
              f"samples {r['tvd_topK_samples_vs_exact_p']['1000']:.4f} (floor {rec['exact_chi']['tvd_topK_sample_floor']['1000']:.4f})  "
              f"unique {r['unique']['unique']} (exact {rec['exact_chi']['unique']['1']['unique']})  "
              f"t_state {info['t_state']:.0f}s sample {t3_ - t2_:.2f}s  iso {smp.isometry_err:.1e}  "
              f"dlogp {r['logp_selfconsistency_max']:.1e}", flush=True)
        del mps, smp
        torch.cuda.empty_cache()
        RES.mkdir(parents=True, exist_ok=True)
        json.dump(rec, open(out_json, "w"), indent=1)              # checkpoint after every chi
        (LOGS / "samples").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_npz, **save)
    del P
    torch.cuda.empty_cache()
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--cand", default="rl4n29f")
    ap.add_argument("--tasks", default="runs_ot/energy_tasks/smallval_rl4n29.pkl")
    ap.add_argument("--chis", type=int, nargs="+", default=[64, 128, 256, 512])
    ap.add_argument("--shots", type=int, default=S.LIN_SHOTS)
    ap.add_argument("--block2-threads", type=int, default=4)
    ap.add_argument("--max-mem-gb", type=float, default=4.0)
    a = ap.parse_args()
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "4")))
    for name in a.names:
        t0 = time.time()
        run_molecule(a, name)
        print(f"[{name}] done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
