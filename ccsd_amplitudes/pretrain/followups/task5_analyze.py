#!/usr/bin/env python3
"""Task 5 analysis: accuracy / cost tables, scaling fits, contention and the NERSC cost-model check (dm TN engine).

Reads pretrain/opt_true/results/followups/task5_dm_scaling/bench.jsonl (task5_dm_bench.py; labels pass1_* = 4
concurrent processes on one GPU, pass2_* = one process alone) and the nvidia-smi logs
rl_runs/followups/task5_dm_scaling/gpu4_monitor_*.jsonl (gpu_monitor.py); writes summary.json next to bench.jsonl and
prints the markdown tables of docs/followups/task5_dm_scaling.md.

    python3 pretrain/followups/task5_analyze.py
"""
from __future__ import annotations

import glob
import importlib.util
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain/opt_true/results/followups/task5_dm_scaling"
LOGD = ROOT / "rl_runs/followups/task5_dm_scaling"
CHIS = (64, 128, 256, 512)
K_CAL = (0.9, 2.2, 6.2)   # err(E) / discarded_sum of the dm engine vs exact, norb 15-18 (tn_energy.py docstring)
PHASES = ("t_state", "t_env", "t_sweep", "t_expect", "t_convert")
BUSY_MAX = 0.25           # a build counts as quiet if the whole GPU was >= 95 % busy for at most this fraction of it


def load_records():
    recs = [json.loads(ln) for ln in open(RES / "bench.jsonl")]
    ok, bad = [r for r in recs if "error" not in r], [r for r in recs if "error" in r]
    for r in ok:
        r["pass"] = r["label"].split("_")[0]
        r["s0"], r["s1"] = r["t_start"], r["t_start"] + r["t_state"]
    for r in ok:            # average number of OTHER benchmark processes building a state during this state build
        r["overlap"] = sum(max(0.0, min(r["s1"], o["s1"]) - max(r["s0"], o["s0"])) for o in ok if o is not r) \
            / max(r["t_state"], 1e-9)
    return ok, bad


def load_monitor():
    rows = []
    for f in sorted(glob.glob(str(LOGD / "gpu4_monitor_*.jsonl"))):
        for ln in open(f):
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if "procs" in r:
                rows.append(r)
    rows.sort(key=lambda r: r["t"])
    return rows


def attach_footprint(recs, mon):
    ts = np.array([r["t"] for r in mon])
    for r in recs:
        i0, i1 = np.searchsorted(ts, r["t_start"]), np.searchsorted(ts, r["t_end"])
        win = mon[i0:i1]
        mem = [w["fts"].get(str(r["pid"]), 0) for w in win]
        r["smi_peak_MiB"] = max(mem) if mem else None
        r["gpu_util_mean"] = float(np.mean([w["gpu_util"] for w in win])) if win else None
        r["foreign_MiB"] = (max(w["gpu_used_MiB"] - sum(w["fts"].values()) for w in win) if win else None)
        r["fts_total_max_MiB"] = max(sum(w["fts"].values()) for w in win) if win else None
        j0, j1 = np.searchsorted(ts, r["s0"]), np.searchsorted(ts, r["s1"])
        u = [w["gpu_util"] for w in mon[j0:j1]]
        # proxy for the foreign job (another user's process shares this GPU): fraction of the build with the whole
        # GPU >= 95 % busy (one dm worker alone keeps it ~40-60 % busy)
        r["busy_frac"] = float(np.mean(np.array(u) >= 95)) if u else None
        r["util_build_mean"] = float(np.mean(u)) if u else None


def fit_power(xs, ys):
    """log y = a log x + b -> (a, stderr(a), n)."""
    x, y = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    if len(set(np.round(x, 9))) < 2:
        return None, None, len(xs)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    se = float("nan")
    if len(x) > 2:
        resid = y - A @ coef
        se = math.sqrt(resid @ resid / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
    return float(coef[0]), se, len(xs)


def med(v):
    v = [x for x in v if x is not None]
    return float(np.median(v)) if v else None


def load_cost_model():
    spec = importlib.util.spec_from_file_location("cm", ROOT / "docs/nersc_2027_cost_model.py")
    cm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cm)
    return cm


def mpo_times(recs):
    """first job of each molecule in each process carries the MPO build."""
    return [r for r in recs if r.get("t_mpo_build", 0.0) > 1.0]


def timing_table(recs):
    tab = defaultdict(lambda: defaultdict(list))
    for r in recs:
        t = tab[(r["norb"], r["chi"])]
        for k in PHASES + ("wall", "smi_peak_MiB", "torch_peak_alloc_MB", "torch_peak_reserved_MB", "overlap",
                           "gpu_util_mean"):
            t[k].append(r.get(k))
        t["names"].append(r["name"])
    mp = defaultdict(list)
    for r in mpo_times(recs):
        mp[r["norb"]].append(r["t_mpo_build"])
    out = {}
    for (n, c), d in sorted(tab.items()):
        row = {k: med(v) for k, v in d.items() if k != "names"}
        row.update({"n_mol": len(d["names"]), "t_mpo_build": med(mp.get(n, [])),
                    "smi_peak_max_MiB": max([x for x in d["smi_peak_MiB"] if x is not None], default=None)})
        out[f"{n}_{c}"] = row
    return out


def fits_vs_norb(recs):
    out = {}
    for chi in CHIS:
        pts = [r for r in recs if r["chi"] == chi]
        for k in ("t_state", "t_env", "t_sweep", "t_expect"):
            if len({r["norb"] for r in pts}) >= 2:
                out[f"{k}_chi{chi}"] = fit_power([r["norb"] / 29 for r in pts], [r[k] for r in pts])
    mp = mpo_times(recs)
    if len({r["norb"] for r in mp}) >= 2:
        out["t_mpo_build"] = fit_power([r["norb"] / 29 for r in mp], [r["t_mpo_build"] for r in mp])
    return out


def chi_exponents(recs):
    per = defaultdict(dict)
    for r in recs:
        per[(r["name"], r["cand"])][r["chi"]] = r
    out = defaultdict(list)
    for (name, _), d in per.items():
        for a, b in ((64, 128), (128, 256), (256, 512)):
            if a in d and b in d:
                for k in ("t_state", "t_expect"):
                    out[f"{k} {a}->{b} n{d[a]['norb']}"].append(math.log(d[b][k] / d[a][k]) / math.log(b / a))
    return {k: (med(v), min(v), max(v), len(v)) for k, v in sorted(out.items())}


def phase_model(recs, chi=128):
    """Per-phase power laws t(n) = t29 * (n/29)^a from the given records at one chi (t29 from the fit)."""
    pts = [r for r in recs if r["chi"] == chi]
    pm = {}
    for k in ("t_state", "t_expect", "t_convert"):
        a, se, m = fit_power([r["norb"] / 29 for r in pts], [r[k] for r in pts])
        x, y = np.log([r["norb"] / 29 for r in pts]), np.log([r[k] for r in pts])
        b = float(np.mean(y - a * x))
        pm[k] = (math.exp(b), a)
    mp = mpo_times(recs)
    a, se, m = fit_power([r["norb"] / 29 for r in mp], [r["t_mpo_build"] for r in mp])
    x, y = np.log([r["norb"] / 29 for r in mp]), np.log([r["t_mpo_build"] for r in mp])
    pm["t_mpo_build"] = (math.exp(float(np.mean(y - a * x))), a)
    return pm


def unit_cost(pm, n, workers=4, per_mpo=8, gpu_bound=False):
    """seconds of one GPU per energy with `workers` workers sharing it: every phase parallel over workers, or the
    MPS build serialized on the GPU (time-slicing) and only the CPU phases overlapping."""
    t = {k: c * (n / 29) ** a for k, (c, a) in pm.items()}
    cpu = t["t_expect"] + t["t_convert"] + t["t_mpo_build"] / per_mpo
    if gpu_bound:
        return max(t["t_state"], (t["t_state"] + cpu) / workers)
    return (t["t_state"] + cpu) / workers


def cost_model_check(cm, pm_scaling):
    """RL lines with the model's n29 A100 anchor but the measured growth with n (ratio unit(n)/unit(29))."""
    c = cm.scenario(1)
    L0, _ = cm.budget(c)
    keys = ["rl_B1", "rl_B2", "rl_B3", "rl_B4", "rl_heavyhex", "rl_dev"]
    base = {k: L0[k] / cm.A100_PER_NODE for k in keys}
    base["rl_total"] = sum(base.values())
    base["total"] = sum(L0.values()) / cm.A100_PER_NODE
    orig = cm.energy_a100_s
    u29 = orig(29, 128, c, cm.GROUP)
    out = {"model_central": base, "model_unit_chi128": {n: orig(n, 128, c, cm.GROUP) for n in (20, 25, 29, 33, 37, 44)}}
    for tag, kw in pm_scaling.items():
        pm, opts = kw

        def patched(n, chi, p, energies_per_mpo, exact_max=19, _pm=pm, _o=opts):
            if n <= exact_max:
                return orig(n, chi, p, energies_per_mpo, exact_max)
            ratio = unit_cost(_pm, n, **_o) / unit_cost(_pm, 29, **_o)
            return orig(29, chi, p, energies_per_mpo, exact_max) * ratio
        cm.energy_a100_s = patched
        try:
            L, _ = cm.budget(c)
        finally:
            cm.energy_a100_s = orig
        row = {k: L[k] / cm.A100_PER_NODE for k in keys}
        row["rl_total"] = sum(row.values())
        row["total"] = sum(L.values()) / cm.A100_PER_NODE
        row["unit_ratio"] = {n: unit_cost(pm, n, **opts) / unit_cost(pm, 29, **opts) for n in (20, 25, 33, 37, 44)}
        row["model_ratio"] = {n: orig(n, 128, c, cm.GROUP) / u29 for n in (20, 25, 33, 37, 44)}
        out[tag] = row
    return out


def main():
    ok, bad = load_records()
    mon = load_monitor()
    attach_footprint(ok, mon)
    S = {"n_ok": len(ok), "n_failed": len(bad),
         "failed": [{k: r.get(k) for k in ("label", "name", "chi", "error")} for r in bad]}

    # ---------------------------------------------------------------- accuracy (E is deterministic: pool the passes)
    acc = {}
    for r in sorted(ok, key=lambda r: r["t_start"]):
        a = acc.setdefault((r["name"], r["cand"]), {"norb": r["norb"], "e_hf": r["e_hf"], "e_ccsd": r["e_ccsd"],
                                                    "E": {}, "disc": {}, "corr": {}, "dmax": {}})
        if r["chi"] in a["E"] and abs(a["E"][r["chi"]] - r["E"]) > 1e-6:
            S.setdefault("reproducibility_warnings", []).append((r["name"], r["chi"], a["E"][r["chi"]], r["E"]))
        a["E"].setdefault(r["chi"], r["E"])
        a["disc"][r["chi"]], a["corr"][r["chi"]], a["dmax"][r["chi"]] = (r["discarded_sum"], r["corr_pct"],
                                                                       r["discarded_max"])
    rows = []
    for (name, cand), a in sorted(acc.items(), key=lambda kv: (kv[1]["norb"], kv[0][1], kv[0][0])):
        chis = sorted(a["E"])
        top = chis[-1]
        ec = (a["e_ccsd"] - a["e_hf"])
        row = {"name": name, "cand": cand, "norb": a["norb"], "chis": chis, "E": a["E"],
               "corr_pct": a["corr"], "disc": a["disc"], "disc_max": a["dmax"],
               "dE_vs_top_mHa": {c: (a["E"][c] - a["E"][top]) * 1e3 for c in chis}, "top": top,
               "ecorr_ccsd_mHa": ec * 1e3,
               "err_est_cal_mHa": {c: [k * a["disc"][c] * 1e3 for k in K_CAL] for c in chis}}
        if len(chis) >= 2:            # E linear in the discarded weight through the two largest chi -> d = 0
            c1, c2 = chis[-2], chis[-1]
            d1, d2, e1, e2 = a["disc"][c1], a["disc"][c2], a["E"][c1], a["E"][c2]
            k = (e1 - e2) / (d1 - d2)
            e0 = e2 - k * d2
            row.update({"extrap_slope_Ha": k, "E_extrap": e0, "err_top_vs_extrap_mHa": (e2 - e0) * 1e3,
                        "err_prev_vs_extrap_mHa": (e1 - e0) * 1e3,
                        "err_top_vs_extrap_pp": (e2 - e0) / -ec * 100})
        rows.append(row)
    S["accuracy"] = rows

    # ---------------------------------------------------------------- timing
    rl = [r for r in ok if r["cand"] == "rl4n29f"]
    best = {}
    for r in rl:                      # fastest build without concurrent own builds; prefer a quiet GPU (busy <= 0.25)
        if r["overlap"] > 0.5:
            continue
        k = (r["name"], r["chi"])
        quiet = r.get("busy_frac") is not None and r["busy_frac"] <= BUSY_MAX
        b = best.get(k)
        bq = b is not None and b.get("busy_frac") is not None and b["busy_frac"] <= BUSY_MAX
        if b is None or (quiet and not bq) or (quiet == bq and r["t_state"] < b["t_state"]):
            best[k] = r
    for r in best.values():
        r["quiet"] = r.get("busy_frac") is not None and r["busy_frac"] <= BUSY_MAX
    cpu = defaultdict(list)
    for r in rl:
        cpu[(r["name"], r["chi"])].append(r)
    mpo_by_mol = defaultdict(list)
    for r in mpo_times(rl):
        mpo_by_mol[r["name"]].append(r["t_mpo_build"])
    clean = []
    for k, r in sorted(best.items(), key=lambda kv: (kv[1]["norb"], kv[0][1], kv[0][0])):
        c = dict(r)
        c["t_expect"] = med([x["t_expect"] for x in cpu[k]])        # CPU phase: not affected by GPU sharing
        c["t_convert"] = med([x["t_convert"] for x in cpu[k]])
        c["t_mpo_build"] = med(mpo_by_mol[r["name"]]) if c["t_mpo_build"] > 1.0 or mpo_by_mol[r["name"]] else 0.0
        c["n_runs"] = len(cpu[k])
        c["t_state_all"] = sorted(round(x["t_state"], 1) for x in cpu[k])
        clean.append(c)
    S["clean_per_job"] = [{kk: c.get(kk) for kk in ("name", "norb", "chi", "t_state", "t_env", "t_sweep", "t_expect",
                                                    "t_convert", "t_mpo_build", "busy_frac", "util_build_mean",
                                                    "overlap", "label", "n_runs", "t_state_all", "smi_peak_MiB", "quiet",
                                                    "torch_peak_alloc_MB", "torch_peak_reserved_MB")}
                          for c in clean]
    quiet = [c for c in clean if c["quiet"]]
    S["timing_clean"] = timing_table(quiet)
    S["fits_clean"] = fits_vs_norb(quiet)
    S["chi_exp_clean"] = chi_exponents(quiet)
    S["timing_clean_incl_busy"] = timing_table(clean)
    for p in ("pass1", "pass2"):
        rr = [r for r in rl if r["pass"] == p]
        if rr:
            S[f"timing_{p}"] = timing_table(rr)
            S[f"fits_{p}"] = fits_vs_norb(rr)
            S[f"chi_exp_{p}"] = chi_exponents(rr)
    # contention: same (molecule, chi) in pass 1 (4 processes) and pass 2 (alone)
    p1 = {(r["name"], r["cand"], r["chi"]): r for r in ok if r["pass"] == "pass1"}
    p2 = {(c["name"], c["cand"], c["chi"]): c for c in clean if c["quiet"]}
    S["contention"] = [{"name": k[0], "cand": k[1], "chi": k[2], "norb": p1[k]["norb"], "overlap_pass1": p1[k]["overlap"],
                        **{f"{ph}_ratio": p1[k][ph] / p2[k][ph] for ph in ("t_state", "t_env", "t_sweep", "t_expect")
                           if p2[k][ph] > 0}}
                       for k in sorted(set(p1) & set(p2), key=lambda k: (p1[k]["norb"], k[2]))]
    # TITAN Xp vs this GPU for identical inputs (C3H5N3_rxn2003_P 'net', tn_v2 n29_dm.jsonl / timing_single.jsonl)
    titan = {}
    for f in ("n29_dm.jsonl", "timing_single.jsonl"):
        for ln in open(ROOT / "pretrain/rl/tests/results/tn_v2" / f):
            r = json.loads(ln)
            if r["name"] == "C3H5N3_rxn2003_P" and r["cand"] == "net":
                titan.setdefault(r["chi"], []).append(r)
    S["titan_vs_here"] = []
    for r in ok:
        if r["cand"] == "net" and r["chi"] in titan:
            t = titan[r["chi"]]
            S["titan_vs_here"].append({"chi": r["chi"], "pass": r["pass"], "E_here": r["E"],
                                       "E_titan": t[0]["E"], "dE_mHa": (r["E"] - t[0]["E"]) * 1e3,
                                       "t_state_here": r["t_state"], "t_state_titan": med([x["t_state"] for x in t]),
                                       "t_expect_here_4thr": r["t_expect"],
                                       "t_expect_titan_1thr": med([x["t_expect"] for x in t]),
                                       "overlap_here": r["overlap"]})

    # ---------------------------------------------------------------- cost model check
    try:
        cm = load_cost_model()
        pm = phase_model(quiet, 128)
        S["phase_model_chi128"] = {k: {"t29_s": v[0], "exponent": v[1]} for k, v in pm.items()}
        S["cost_model"] = cost_model_check(cm, {
            "measured_scaling_all_parallel": (pm, {"workers": 4, "per_mpo": 8}),
            "measured_scaling_gpu_bound_build": (pm, {"workers": 4, "per_mpo": 8, "gpu_bound": True}),
            "measured_scaling_mpo_every_task": (pm, {"workers": 4, "per_mpo": 1}),
            "measured_scaling_gpu_bound_mpo_every_task": (pm, {"workers": 4, "per_mpo": 1, "gpu_bound": True}),
        })
        # RL line for an n29 anchor A (A100-s per chi-128 energy with 4 workers) and the measured GPU-bound growth
        u29_model = S["cost_model"]["model_unit_chi128"][29]
        tn_keys = ("rl_B2", "rl_B3", "rl_B4", "rl_heavyhex", "rl_dev")
        b1 = S["cost_model"]["model_central"]["rl_B1"]
        S["rl_line_vs_anchor"] = {}
        for growth in ("measured_scaling_gpu_bound_mpo_every_task", "measured_scaling_all_parallel",
                       "measured_scaling_mpo_every_task"):
            gb = S["cost_model"][growth]
            tn_line = sum(gb[k] for k in tn_keys)
            S["rl_line_vs_anchor"][growth] = {f"{a:.1f}": b1 + tn_line * a / u29_model
                                              for a in (28.7, 34.8, 40.0, 48.8)}
        S["rl_line_vs_anchor_note"] = ("node-hours of the RL lines (B1 exact unchanged); anchor = A100-s per n29 chi-128 "
                                       "energy with 4 workers; 28.7 / 48.8 = measured with 4 workers on the shared H200 / a TITAN Xp (task 4); 34.8 = model")
        S["unit_cost_here_s"] = {tag: {n: unit_cost(pm, n, **o) for n in (29, 33, 37, 44)} for tag, o in (
            ("all_parallel_mpo8", {"workers": 4, "per_mpo": 8}), ("gpu_bound_mpo8", {"workers": 4, "per_mpo": 8, "gpu_bound": True}),
            ("gpu_bound_mpo_every_task", {"workers": 4, "per_mpo": 1, "gpu_bound": True}),
            ("single_worker_mpo_every_task", {"workers": 1, "per_mpo": 1}))}
    except Exception as e:  # noqa: BLE001
        import traceback
        S["cost_model_error"] = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"

    # ---------------------------------------------------------------- 4-worker queue throughput (RL configuration)
    qt = []
    f = ROOT / "pretrain/opt_true/results/bench_dm_chi128_h200_4x4.json"          # task 4 (coordinator), n29
    if f.exists():
        d = json.load(open(f))
        nt = sum(len(v) for v in d["per_molecule"].values())
        qt.append({"norb": 29, "n_tasks": nt, "wall_s": d["wall_s"], "s_per_energy": d["wall_s"] / nt,
                   "energies_per_h": 3600 * nt / d["wall_s"], "source": "task 4 (bench_dm_chi128_h200_4x4.json)"})
    f = RES / "queue_tput.json"
    if f.exists():
        for b in json.load(open(f))["batches"]:
            qt.append({k: b[k] for k in ("norb", "n_tasks", "wall_s", "s_per_energy", "energies_per_h")}
                      | {"source": "task 5 (queue_tput.json)", "task_t_median": b["task_t_median"]})
    S["queue_throughput"] = qt
    json.dump(S, open(RES / "summary.json", "w"), indent=1, default=str)
    print_tables(S)
    if qt:
        print("\n### 4 workers x 4 threads through the queue (dm, chi 128; MPO rebuilt for most tasks)\n")
        base = next((q for q in qt if q["norb"] == 29), None)
        print("| norb | tasks | wall (s) | s per energy | energies / GPU-h | ratio to n29 | source |")
        print("|--:|--:|--:|--:|--:|--:|:--|")
        for q in qt:
            r = q["s_per_energy"] / base["s_per_energy"] if base else float("nan")
            print(f"| {q['norb']} | {q['n_tasks']} | {q['wall_s']:.0f} | {q['s_per_energy']:.1f} | {q['energies_per_h']:.0f} | "
                  f"{r:.2f} | {q['source']} |")


def print_tables(S):
    print(f"records ok {S['n_ok']}, failed {S['n_failed']}")
    for f in S["failed"]:
        print("  FAILED", f)
    for w in S.get("reproducibility_warnings", []):
        print("  NOT REPRODUCED", w)
    print("\n### accuracy: % of the CCSD correlation energy, ΔE vs the largest χ (mHa), discarded_sum\n")
    print("| molecule | norb | cand | " + " | ".join(f"% χ{c}" for c in CHIS) + " | " +
          " | ".join(f"ΔE χ{c}" for c in CHIS[:-1]) + " | " + " | ".join(f"disc χ{c}" for c in CHIS) +
          " | χtop − extrap (mHa) | slope (Ha) |")
    print("|:--|--:|:--|" + "--:|" * (3 * len(CHIS) - 1 + 2))
    for r in S["accuracy"]:
        f = lambda d, c, fmt: (fmt % d[c]) if c in d else "–"  # noqa: E731
        print(f"| {r['name']} | {r['norb']} | {r['cand']} | " + " | ".join(f(r['corr_pct'], c, '%.2f') for c in CHIS)
              + " | " + " | ".join(f(r['dE_vs_top_mHa'], c, '%.2f') for c in CHIS[:-1]) + " | "
              + " | ".join(f(r['disc'], c, '%.1e') for c in CHIS)
              + f" | {r.get('err_top_vs_extrap_mHa', float('nan')):.2f} | {r.get('extrap_slope_Ha', float('nan')):.2f} |")
    for p in ("clean", "pass2", "pass1"):
        tab = S.get(f"timing_{p}")
        if not tab:
            continue
        print(f"\n### timing medians, {p} (s; memory MiB; overlap = mean number of other concurrent state builds)\n")
        print("| norb | χ | mols | t_state | t_env | t_sweep | t_expect (4 thr) | t_convert | t_mpo_build | "
              "nvidia-smi peak (max) | torch alloc | torch reserved | overlap | GPU util % |")
        print("|--:|--:|--:|" + "--:|" * 11)
        for key, d in tab.items():
            n, c = key.split("_")
            g = lambda k, fmt="%.1f": (fmt % d[k]) if d.get(k) is not None else "–"  # noqa: E731
            print(f"| {n} | {c} | {d['n_mol']} | {g('t_state')} | {g('t_env')} | {g('t_sweep')} | {g('t_expect')} | "
                  f"{g('t_convert')} | {g('t_mpo_build')} | {g('smi_peak_MiB', '%.0f')} ({g('smi_peak_max_MiB', '%.0f')}) | "
                  f"{g('torch_peak_alloc_MB', '%.0f')} | {g('torch_peak_reserved_MB', '%.0f')} | "
                  f"{g('overlap', '%.1f')} | {g('gpu_util_mean', '%.0f')} |")
        fits = S.get(f"fits_{p}", {})
        print(f"\nexponents vs norb ({p}; fit of log t vs log n over molecules): " + "; ".join(
            f"{k} {v[0]:.2f} ± {v[1]:.2f} (n={v[2]})" for k, v in fits.items() if v[0] is not None))
        cx = S.get(f"chi_exp_{p}", {})
        print(f"\nχ exponents ({p}, median [min, max] over molecules): " + "; ".join(
            f"{k}: {v[0]:.2f} [{v[1]:.2f}, {v[2]:.2f}]" for k, v in cx.items()))
    if S.get("contention"):
        print("\n### contention: pass-1 time (4 processes) / quiet single-process time, same job\n")
        print("| molecule | norb | χ | overlap in pass 1 | t_state | t_env | t_sweep | t_expect |")
        print("|:--|--:|--:|--:|--:|--:|--:|--:|")
        for c in S["contention"]:
            print(f"| {c['name']} | {c['norb']} | {c['chi']} | {c['overlap_pass1']:.1f} | {c.get('t_state_ratio', 0):.2f} | "
                  f"{c.get('t_env_ratio', 0):.2f} | {c.get('t_sweep_ratio', 0):.2f} | {c.get('t_expect_ratio', 0):.2f} |")
    if S.get("titan_vs_here"):
        print("\n### identical inputs (C3H5N3_rxn2003_P, 'net' parameters): TITAN Xp (tn_v2) vs this run\n")
        for t in S["titan_vs_here"]:
            print(t)
    if "phase_model_chi128" in S:
        print("\nphase model at χ128 (t29, exponent):", {k: (round(v['t29_s'], 1), round(v['exponent'], 2))
                                                        for k, v in S["phase_model_chi128"].items()})
    if "unit_cost_here_s" in S:
        print("\nper-energy seconds of one GPU on this node (4 workers unless stated):",
              json.dumps(S["unit_cost_here_s"], default=lambda x: round(x, 1)))
    if "rl_line_vs_anchor" in S:
        print("\nRL line (node-h) vs the n29 anchor, measured GPU-bound growth:", json.dumps(S["rl_line_vs_anchor"],
              default=lambda x: round(x)), "--", S["rl_line_vs_anchor_note"])
    if "cost_model" in S:
        print("\n### cost model RL lines (node-hours)\n")
        for k, v in S["cost_model"].items():
            print(k, json.dumps(v, default=lambda x: round(x, 2)))
    if "cost_model_error" in S:
        print(S["cost_model_error"])


if __name__ == "__main__":
    sys.exit(main())
