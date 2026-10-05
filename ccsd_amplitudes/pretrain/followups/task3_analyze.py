#!/usr/bin/env python3
"""Task 3 analysis: 500-iteration optimize=True labels vs the same fits run to convergence vs the network candidates.

Reads the existing result JSONs (merged as in pretrain/rl/report_tables.py), the new energies in
pretrain/opt_true/results/followups/task3_converged_labels/, and the fit metadata in rhf_targets_compressed_conv/.
Prints markdown tables (stdout) and writes summary.json next to the new energies.

Candidates
    label     the existing label (L-BFGS capped at 500 iterations) and its stored energy
    lab500    this run's own 500-iteration point (snapshot), on the trajectory that was continued to convergence.
              Identical to `label` where the trajectory reproduces the stored label bit-for-bit; evaluated anew where
              it does not (norb 19 and the n29 test set: those labels were made in a run whose floating-point
              reduction order differed, and the L-BFGS trajectory is chaotic)
    labconv   the same fit run to convergence (scipy L-BFGS-B success, maxiter 5000)

    pretrain/.train_venv/bin/python3 pretrain/followups/task3_analyze.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain" / "opt_true" / "results"
NEW = RES / "followups" / "task3_converged_labels"
CONV = ROOT / "rhf_targets_compressed_conv" / "square_reg0.005"

SETS = {  # title: (names file, existing result files, label dirs behind the stored label energies, new files, kind)
    "norb 15-16 val (9)": ("pretrain/rl/small_val.txt",
                           ["energy_small_baselines.json", "energy_small_all4.json", "energy_smallval_rl4L.json",
                            "energy_smallval_rl4n29.json"],
                           ["rhf_targets_compressed_small"], ["energy_smallval_labconv.json"], "exact"),
    "norb 17-18 val (22)": ("pretrain/rl/large_val.txt",
                            ["energy_largeval_base.json", "energy_largeval_pre4_rl4L.json",
                             "energy_largeval_rl4n29.json", "energy_largeval_rl4n29_C2H6O.json"],
                                  ["rhf_targets_compressed_n17", "rhf_targets_compressed_n18"],
                            ["energy_largeval_labconv.json", "energy_largeval_lab500.json"], "exact"),
    "norb 19 test (12)": ("gauge_study/names_norb19.txt", ["energy_n19_all.json", "energy_n19_rl4n29.json"],
                          ["rhf_targets_compressed_n19"], ["energy_n19_labconv.json", "energy_n19_warm.json", "energy_n19_warmown.json"],
                          "exact"),
    "n29 val (20)": ("pretrain/rl/n29_val20.txt",
                     ["energy_n29val_base_chi256.json", "energy_n29val_rl4n29_chi256.json"],
                     ["rhf_targets_compressed"], ["energy_n29val_labconv_chi256.json"], "MPS chi 256, v1 engine"),
    "n29 test (24)": ("pretrain/rl/n29_test24.txt",
                      ["energy_n29test_base_chi256.json", "energy_n29test_rl4n29_chi256.json"],
                      ["rhf_targets_compressed_n29test"],
                      ["energy_n29test_labconv_chi256.json", "energy_n29test_lab500_chi256.json",
                       "energy_n29test_warm_chi256.json"],
                      "MPS chi 256, v1 engine"),
}
RENAME = {"all4_t4": "pre4"}          # energy_small_all4.json: the pretrained slot model at 4 recycles
NET = ["pre4", "rl4L", "rl4n29", "rl4n29f"]
LAB = {"label": "label, 500 it. (stored)", "lab500": "label, 500 it. (this trajectory)",
       "labconv": "label, converged", "warm_conv": "stored label continued to convergence (warm start)",
       "warmown_conv": "this trajectory's 500-it. label continued to convergence (warm start)", "pre4": "pretrained, 4 recycles", "rl4L": "norb 15-18 RL (step 25)",
       "rl4n29": "+ n29 TN RL (step 15)", "rl4n29f": "+ n29 TN RL (step 30)"}
TIE = 5e-7     # |difference| < 5e-5 points of % corr (~1e-7 Ha): a tie. Identical parameters re-evaluated by the
#                complex64 exact engine differ by <= 1.5e-7 here; different states never came closer than 1.5e-6


def load_pm(path, keep_E=False):
    p = Path(path)
    if not p.exists():
        return {}, {}
    d = json.load(open(p))
    pm = {}
    for n, cs in d["per_molecule"].items():
        for c, v in cs.items():
            if v.get("error") or v.get("E") is None or not np.isfinite(v["E"]):
                continue
            pm.setdefault(n, {})[RENAME.get(c, c)] = v if keep_E else v["corr_frac"]
    meta = {n: {RENAME.get(c, c): r for c, r in m.get("resid", {}).items()} for n, m in d.get("meta", {}).items()}
    return pm, meta


def merge(files, base, keep_E=False):
    per, res = {}, {}
    for f in files:
        pm, me = load_pm(base / f, keep_E)
        for n, d in pm.items():
            per.setdefault(n, {}).update(d)
        for n, d in me.items():
            res.setdefault(n, {}).update(d)
    return per, res


def fit_meta(name, label_dirs):
    """Fit summary of the long run + bit-for-bit check against the stored label behind the stored energies."""
    z = np.load(CONV / f"{name}.npz")
    o = next(np.load(p) for p in (ROOT / d / "square_reg0.005" / f"{name}.npz" for d in label_dirs) if p.exists())
    its = [int(k) for k in z["snap_iters"]]
    k_old = int(o["nit"])
    if k_old in its:
        j = its.index(k_old)
        Zs, Us = z["snap_Z"][j], z["snap_U_re"][j] + 1j * z["snap_U_im"][j]
    elif k_old == int(z["nit"]):
        Zs, Us = z["Z"], z["U_re"] + 1j * z["U_im"]
    else:
        Zs = Us = None
    Uo = o["U_re"] + 1j * o["U_im"]
    dU = float(np.max(np.abs(Us - Uo))) if Us is not None else float("nan")
    dZ = float(np.max(np.abs(Zs - o["Z"]))) if Zs is not None else float("nan")
    r500 = float(z["snap_resid"][its.index(500)]) if 500 in its else float(z["resid"])
    return {"nit": int(z["nit"]), "nfev": int(z["nfev"]), "success": bool(z["success"]),
            "message": str(z["message"]), "resid": float(z["resid"]), "gmax": float(z["gmax"]),
            "time_s": float(z["time"]), "resid500_own": r500, "old_nit": k_old, "old_success": bool(o["success"]),
            "old_resid": float(o["resid"]), "traj_dU": dU, "traj_dZ": dZ, "same_traj": bool(dU == 0.0 and dZ == 0.0)}


def stats(v):
    v = np.asarray(v)
    return {"n": int(len(v)), "mean": float(100 * v.mean()), "median": float(100 * np.median(v)),
            "min": float(100 * v.min())}


def paired(a, b):
    d = np.asarray(a) - np.asarray(b)
    return {"mean": float(100 * d.mean()), "median": float(100 * np.median(d)), "wins": int((d > TIE).sum()),
            "losses": int((d < -TIE).sum()), "ties": int((np.abs(d) <= TIE).sum()),
            "n": int(len(d)), "min": float(100 * d.min()), "max": float(100 * d.max())}


def wl(p):
    return f"{p['wins']}/{p['n']}" + (f" ({p['ties']} tied)" if p["ties"] else "")


def spearman(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3:
        return float("nan")
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def spread_section(summary):
    """Run-to-run spread of the label at norb 15-17 (exact energies). Runs with identical inputs and settings:
    p48 = the stored labels (multi-threaded XLA; the 2-thread run is bit-identical to it, so it is not counted),
    p1 = single-threaded XLA (different float summation order), s1-s4 = canonical init x0 * (1 + 1e-15 N(0,1))."""
    names = [ln.strip() for ln in open(ROOT / "pretrain/followups/task3_lists/traj_n15_17.txt") if ln.strip()]
    sm, _ = merge(SETS["norb 15-16 val (9)"][1], RES)
    lg, _ = merge(SETS["norb 17-18 val (22)"][1], RES)
    old = {**sm, **lg}
    main, _ = merge(["energy_smallval_labconv.json", "energy_largeval_labconv.json"], NEW)
    runs = ["p1"] + [f"s{k}" for k in range(1, 5)]
    var, _ = merge([f"energy_traj_{t}.json" for t in runs], NEW)
    rows = []
    for n in names:
        e500, econv = {"p48": old[n].get("label")}, {"p48": main.get(n, {}).get("labconv")}
        for t in runs:
            fp = ROOT / "rhf_targets_compressed_conv" / f"traj_{t}" / "square_reg0.005" / f"{n}.npz"
            if not fp.exists():
                continue
            nit = int(np.load(fp)["nit"])
            econv[t] = var.get(n, {}).get(f"{t}_conv")
            e500[t] = econv[t] if nit <= 500 else var.get(n, {}).get(f"{t}_500")
        if any(v is None for v in list(e500.values()) + list(econv.values())):
            continue
        rows.append((n, e500, econv, {c: old[n].get(c) for c in NET}))
    if not rows:
        return
    R = sorted(rows[0][2])
    rows = [r for r in rows if sorted(r[2]) == R]
    K = len(R)
    print(f"## Label run-to-run spread, norb 15-17 ({len(rows)} molecules x {K} runs, exact energies)\n")
    print("Identical inputs and settings (canonical init, square_reg0.005, default L-BFGS-B tolerances, maxiter 5000). "
          "Runs: the stored labels (multi-threaded XLA), a single-threaded run (different floating-point summation "
          "order), and 4 runs from the canonical init perturbed by relative 1e-15 Gaussian noise.\n")
    A = {k: np.array([[r[i][t] for t in R] for r in rows]) for k, i in (("500", 1), ("conv", 2))}   # (mol, run)
    print("| quantity (points of % corr) | 500 iterations | converged |")
    print("|:--|--:|--:|")
    rng = {k: 100 * (A[k].max(1) - A[k].min(1)) for k in A}
    sd = {k: 100 * A[k].std(1, ddof=1) for k in A}
    setm = {k: 100 * A[k].mean(0) for k in A}
    print(f"| per-molecule SD over the {K} runs: mean / median / max | {sd['500'].mean():.2f} / "
          f"{np.median(sd['500']):.2f} / {sd['500'].max():.2f} | {sd['conv'].mean():.2f} / {np.median(sd['conv']):.2f} / "
          f"{sd['conv'].max():.2f} |")
    print(f"| per-molecule range (best − worst run): mean / max | {rng['500'].mean():.2f} / {rng['500'].max():.2f} | "
          f"{rng['conv'].mean():.2f} / {rng['conv'].max():.2f} |")
    print(f"| molecules where all {K} runs agree (range < 0.01) | {int((rng['500'] < 0.01).sum())}/{len(rows)} | "
          f"{int((rng['conv'] < 0.01).sum())}/{len(rows)} |")
    print(f"| set mean per run: stored run / min / max over runs | {setm['500'][R.index('p48')]:.2f} / "
          f"{setm['500'].min():.2f} / {setm['500'].max():.2f} | {setm['conv'][R.index('p48')]:.2f} / "
          f"{setm['conv'].min():.2f} / {setm['conv'].max():.2f} |")
    print(f"| SD of the set mean over runs | {setm['500'].std(ddof=1):.2f} | {setm['conv'].std(ddof=1):.2f} |")
    print(f"| set mean of the per-molecule best run (selected by energy) | {100 * A['500'].max(1).mean():.2f} | "
          f"{100 * A['conv'].max(1).mean():.2f} |")
    print()
    print(f"| network | mean % corr | vs mean converged label | P(network > a converged label) | "
          f"> best of {K} converged labels | vs best of {K} (mean) |")
    print("|:--|--:|--:|--:|--:|--:|")
    net = {}
    for c in NET:
        v = np.array([r[3][c] if r[3][c] is not None else np.nan for r in rows])
        if np.isnan(v).any():
            continue
        C = A["conv"]
        net[c] = {"mean": float(100 * v.mean()), "vs_mean_label": float(100 * (v - C.mean(1)).mean()),
                  "p_beats": float((v[:, None] > C).mean()), "beats_best": int((v > C.max(1)).sum()),
                  "vs_best": float(100 * (v - C.max(1)).mean())}
        print(f"| {LAB[c]} | {net[c]['mean']:.2f} | {net[c]['vs_mean_label']:+.2f} | {net[c]['p_beats']:.2f} | "
              f"{net[c]['beats_best']}/{len(v)} | {net[c]['vs_best']:+.2f} |")
    print()
    summary["label_spread_norb15_17"] = {
        "runs": R, "sd_points": {k: sd[k].tolist() for k in sd}, "range_points": {k: rng[k].tolist() for k in rng},
        "set_means": {k: dict(zip(R, setm[k].tolist())) for k in setm}, "network": net,
        "per_molecule": {r[0]: {"e500": r[1], "econv": r[2]} for r in rows}}


def tight_section(summary):
    """Default L-BFGS-B stop (gtol 1e-5) vs gtol 1e-10 / ftol 1e-15 on the SAME trajectory (both single-threaded)."""
    pairs = []
    sm_t, _ = merge(["energy_smallval_tight.json"], NEW)
    p1, _ = merge(["energy_traj_p1.json"], NEW)
    for n, d in sm_t.items():
        if "tight_conv" in d and "p1_conv" in p1.get(n, {}):
            a = np.load(ROOT / "rhf_targets_compressed_conv/traj_p1/square_reg0.005" / f"{n}.npz")
            pairs.append((n, "exact", a, d["tight_conv"], p1[n]["p1_conv"]))
    nt, _ = merge(["energy_n29val_tight_chi256.json"], NEW)
    nc, _ = merge(["energy_n29val_labconv_chi256.json"], NEW)
    for n, d in nt.items():
        if "tight_conv" in d and "labconv" in nc.get(n, {}):
            a = np.load(CONV / f"{n}.npz")
            pairs.append((n, "MPS chi 256", a, d["tight_conv"], nc[n]["labconv"]))
    if not pairs:
        return
    print(f"## Tighter L-BFGS tolerance on the same trajectory ({len(pairs)} molecules)\n")
    print("| molecule | energies | iterations default → tight | t2 resid default → tight | % corr default → tight |")
    print("|:--|:--|--:|--:|--:|")
    rows = []
    for n, kind, a, et, ed in pairs:
        b = np.load(ROOT / "rhf_targets_compressed_conv/tight/square_reg0.005" / f"{n}.npz")
        rows.append({"name": n, "kind": kind, "nit": [int(a["nit"]), int(b["nit"])],
                     "resid": [float(a["resid"]), float(b["resid"])], "corr": [100 * ed, 100 * et]})
        print(f"| {n} | {kind} | {int(a['nit'])} → {int(b['nit'])} | {float(a['resid']):.4f} → {float(b['resid']):.4f} "
              f"| {100 * ed:.2f} → {100 * et:.2f} ({100 * (et - ed):+.2f}) |")
    d = np.array([r["corr"][1] - r["corr"][0] for r in rows])
    print(f"\nEnergy change from tightening: mean {d.mean():+.2f}, mean |Δ| {np.abs(d).mean():.2f}, "
          f"max |Δ| {np.abs(d).max():.2f} points.\n")
    summary["tight_vs_default"] = rows


def main():
    summary, pooled = {}, []
    for title, (nf, old_files, label_dirs, new_files, kind) in SETS.items():
        names = [ln.strip() for ln in open(ROOT / nf) if ln.strip()]
        per, res = merge(old_files, RES)
        new, new_res = merge(new_files, NEW)
        newE, _ = merge(new_files, NEW, keep_E=True)
        oldE, _ = merge(old_files, RES, keep_E=True)
        fm = {n: fit_meta(n, label_dirs) for n in names if (CONV / f"{n}.npz").exists()}
        for n in names:
            per.setdefault(n, {})
            res.setdefault(n, {})
            for c in ("labconv", "lab500", "label_ctl", "warm_conv", "warmown_conv"):
                if c in new.get(n, {}):
                    per[n][c] = new[n][c]
                if c in new_res.get(n, {}):
                    res[n][c] = new_res[n][c]
            if n in fm:
                res[n]["lab500"] = fm[n]["resid500_own"]
                res[n]["labconv"] = fm[n]["resid"]
                if "lab500" not in per[n]:
                    if fm[n]["nit"] <= 500 and "labconv" in per[n]:
                        per[n]["lab500"] = per[n]["labconv"]       # converged before 500: same parameters
                    elif fm[n]["same_traj"] and fm[n]["old_nit"] == 500 and "label" in per[n]:
                        per[n]["lab500"] = per[n]["label"]         # bit-identical to the stored label
        ok = [n for n in names if all(c in per[n] for c in ("label", "lab500", "labconv"))]
        print(f"## {title}: {kind} energies, {len(ok)}/{len(names)} molecules complete\n")
        if not ok:
            print("(pending)\n")
            continue
        f = [fm[n] for n in ok]
        nit = np.array([x["nit"] for x in f])
        same = sum(x["same_traj"] for x in f)
        print(f"* Fits (maxiter 5000): {sum(x['success'] for x in f)}/{len(f)} converged "
              f"({'; '.join(sorted({x['message'] for x in f}))}); iterations median {int(np.median(nit))}, "
              f"range {nit.min()}-{nit.max()}; {np.median([x['time_s'] for x in f]):.0f} s per fit (median, 1 core). "
              f"Stored labels: {sum(x['old_success'] for x in f)}/{len(f)} converged, "
              f"{sum(x['old_nit'] >= 500 for x in f)} stopped at the 500 cap.")
        print(f"* Same trajectory as the stored labels (bit-identical parameters at the stored label's iteration): "
              f"{same}/{len(f)}" + ("" if same == len(f) else
                                     f"; the others differ by max |ΔU| {np.nanmedian([x['traj_dU'] for x in f if not x['same_traj']]):.2f} "
                                     f"(median), so `lab500` was evaluated on this run's trajectory") + ".\n")
        cands = (["label"] + (["lab500"] if same < len(f) else []) + ["labconv"]
                 + [c for c in ("warm_conv", "warmown_conv") if any(c in per[n] for n in ok)] + NET)
        lab = np.array([per[n]["label"] for n in ok])
        conv = np.array([per[n]["labconv"] for n in ok])
        print("| candidate | n | mean t2 resid | mean % corr | median | min | vs stored label | > stored label "
              "| vs converged label | > converged label |")
        print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
        tab = {}
        for c in cands:
            idx = [i for i, n in enumerate(ok) if c in per[n]]
            if not idx:
                continue
            v = np.array([per[ok[i]][c] for i in idx])
            rr = [res[ok[i]].get(c) for i in idx]
            rmean = f"{np.mean(rr):.3f}" if all(r is not None for r in rr) else ""
            s = stats(v)
            p5, pc = paired(v, lab[idx]), paired(v, conv[idx])
            tab[c] = {**s, "resid_mean": None if not rmean else float(rmean), "vs_label": p5, "vs_labconv": pc}
            a, aw = ("—", "—") if c == "label" else (f"{p5['mean']:+.2f}", wl(p5))
            b, bw = ("—", "—") if c == "labconv" else (f"{pc['mean']:+.2f}", wl(pc))
            print(f"| {LAB.get(c, c)} | {s['n']} | {rmean} | {s['mean']:.2f} | {s['median']:.2f} | {s['min']:.2f} | "
                  f"{a} | {aw} | {b} | {bw} |")
        print()
        for alt in ("warm_conv", "warmown_conv"):
            if all(alt in per[n] for n in ok):
                av = np.array([per[n][alt] for n in ok])
                txt = []
                for c in NET:
                    if all(c in per[n] for n in ok):
                        v = np.array([per[n][c] for n in ok])
                        p_ = paired(v, av)
                        txt.append(f"{c} {p_['mean']:+.2f} ({wl(p_)})")
                print(f"* Networks vs `{alt}`: " + "; ".join(txt) + ".")
                tab.setdefault("vs_" + alt, {c: paired(np.array([per[n][c] for n in ok]), av) for c in NET
                                             if all(c in per[n] for n in ok)})
        alts = ["labconv"] + [c for c in ("warm_conv", "warmown_conv") if all(c in per[n] for n in ok)]
        if len(alts) > 1:
            best = np.max(np.array([[per[n][c] for c in alts] for n in ok]), 1)
            txt = [f"{c} {paired(np.array([per[n][c] for n in ok]), best)['mean']:+.2f} "
                   f"({wl(paired(np.array([per[n][c] for n in ok]), best))})" for c in NET if all(c in per[n] for n in ok)]
            print(f"* Networks vs the best of the {len(alts)} converged labels per molecule ({', '.join(alts)}; set mean "
                  f"{100 * best.mean():.2f}): " + "; ".join(txt) + ".\n")
            tab["best_conv"] = {"alts": alts, "mean": float(100 * best.mean()),
                                **{c: paired(np.array([per[n][c] for n in ok]), best) for c in NET
                                   if all(c in per[n] for n in ok)}}
        # same trajectory: 500 iterations -> converged
        e500 = np.array([per[n]["lab500"] for n in ok])
        r500 = np.array([fm[n]["resid500_own"] for n in ok])
        rc = np.array([fm[n]["resid"] for n in ok])
        de, dr = 100 * (conv - e500), rc - r500
        mv = nit > 500
        print(f"* 500 iterations -> converged, same trajectory ({int(mv.sum())}/{len(ok)} fits ran past 500): "
              f"t2 residual change mean {dr[mv].mean() if mv.any() else 0:+.5f}, median "
              f"{np.median(dr[mv]) if mv.any() else 0:+.5f} (range {dr.min():+.4f} .. {dr.max():+.4f}); energy change mean {de[mv].mean() if mv.any() else 0:+.2f}, median "
              f"{np.median(de[mv]) if mv.any() else 0:+.2f} points of % corr (range {de.min():+.2f} .. {de.max():+.2f}); "
              f"better on {int((de[mv] > 100 * TIE).sum())}, worse on {int((de[mv] < -100 * TIE).sum())}; "
              f"Spearman(Δresid, Δ% corr) {spearman(dr[mv], de[mv]):+.2f}.")
        # trajectory noise: stored label vs this run's 500-iteration point (both 500 iterations, same settings)
        dif = [n for n in ok if not fm[n]["same_traj"]]
        if dif:
            dn = np.array([100 * (per[n]["label"] - per[n]["lab500"]) for n in dif])
            drn = np.array([fm[n]["old_resid"] - fm[n]["resid500_own"] for n in dif])
            dcv = np.array([100 * (per[n]["label"] - per[n]["labconv"]) for n in dif])
            print(f"* Trajectory noise ({len(dif)} molecules whose stored label came from a different floating-point "
                  f"trajectory): stored label − this run at 500 iterations = {dn.mean():+.2f} points mean, "
                  f"mean |Δ| {np.abs(dn).mean():.2f}, max |Δ| {np.abs(dn).max():.2f}; residual difference mean |Δ| "
                  f"{np.abs(drn).mean():.4f}; stored label − converged label mean |Δ| {np.abs(dcv).mean():.2f}.")
        ctl = []
        for n in ok:
            if "label_ctl" in newE.get(n, {}) and "label" in oldE.get(n, {}):
                ctl.append((n, 1e3 * (newE[n]["label_ctl"]["E"] - oldE[n]["label"]["E"])))
            elif fm[n]["same_traj"] and fm[n]["nit"] == fm[n]["old_nit"] and "labconv" in newE.get(n, {}):
                ctl.append((n, 1e3 * (newE[n]["labconv"]["E"] - oldE[n]["label"]["E"])))
        if ctl:
            print(f"* Controls (identical parameters re-evaluated now vs the stored energy): {len(ctl)} molecule(s), "
                  f"max |ΔE| {max(abs(x[1]) for x in ctl):.4f} mHa.")
        print()
        summary[title] = {"kind": kind, "n": len(ok), "table": tab, "fits": {n: fm[n] for n in ok},
                          "per_molecule": {n: {c: per[n][c] for c in cands + ["lab500", "label_ctl"] if c in per[n]}
                                           for n in ok},
                          "resid": {n: {c: res[n][c] for c in res[n]} for n in ok},
                          "same_traj_500_to_conv": {"dresid": dr.tolist(), "dcorr_points": de.tolist()},
                          "controls_mHa": ctl}
        pooled += [(title, n, per[n]) for n in ok]
    if pooled:
        print(f"## All sets pooled ({len(pooled)} molecules; exact energies at norb <= 19, MPS chi 256 at n29)\n")
        print("| candidate | n | mean % corr | vs stored label | > stored label | vs converged label "
              "| > converged label |")
        print("|:--|--:|--:|--:|--:|--:|--:|")
        out = {}
        for c in ["label", "lab500", "labconv"] + NET:
            rows = [(p["label"], p["labconv"], p[c]) for (_, _, p) in pooled if c in p]
            lb, cv, x = map(np.array, zip(*rows))
            p5, pc = paired(x, lb), paired(x, cv)
            out[c] = {"n": len(x), "mean": float(100 * x.mean()), "vs_label": p5, "vs_labconv": pc}
            print(f"| {LAB.get(c, c)} | {len(x)} | {100 * x.mean():.2f} | {p5['mean']:+.2f} | {wl(p5)} | "
                  f"{pc['mean']:+.2f} | {wl(pc)} |")
        summary["pooled"] = out
        print()
    spread_section(summary)
    tight_section(summary)
    NEW.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(NEW / "summary.json", "w"), indent=1)
    print(f"-> {NEW.relative_to(ROOT) / 'summary.json'}")


if __name__ == "__main__":
    main()
