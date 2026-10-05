#!/usr/bin/env python3
"""Task 3, section 4 of docs/followups/task3_converged_labels.md: does the t2 residual of a label predict its energy?

Reads summary.json (written by task3_analyze.py) and the fit .npz files in rhf_targets_compressed_conv/; prints the
statistics and writes summary_resid_vs_energy.json next to the energies.

  (a) six equivalent label runs per molecule (norb 15-17, 21 molecules; the runs of task3_analyze.spread_section):
      per-molecule Spearman(resid, % corr), pairwise concordance (lower resid => better energy), set mean when the
      lowest- / highest-residual run is picked, run mean, per-molecule best run (oracle)
  (b) all label variants of the same molecule (label, lab500, labconv, warm_conv, warmown_conv), 87 molecules
  (c) networks vs the converged label: residual higher and energy better
  (d) across molecules, converged labels: Spearman(resid, % corr) per set (molecule difficulty)

    pretrain/.train_venv/bin/python3 pretrain/followups/task3_resid_energy.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
NEW = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "task3_converged_labels"
CONV = ROOT / "rhf_targets_compressed_conv"
LABS = ["label", "lab500", "labconv", "warm_conv", "warmown_conv"]
NET = ["pre4", "rl4L", "rl4n29", "rl4n29f"]


def main():
    d = json.load(open(NEW / "summary.json"))
    sets = [k for k in d if k not in ("pooled", "label_spread_norb15_17", "tight_vs_default")]
    res, ene = {}, {}
    for s in sets:
        for n, v in d[s]["resid"].items():
            res.setdefault(n, {}).update(v)
        for n, v in d[s]["per_molecule"].items():
            ene.setdefault(n, {}).update(v)
    out = {"six_runs": {}, "label_variants": {}, "networks_vs_labconv": {}, "across_molecules_labconv": {}}

    # (a) six equivalent runs
    sp = d["label_spread_norb15_17"]
    runs = sp["runs"]
    for stage in ("500", "conv"):
        rhos, hit, pick_lo, pick_hi, mean, best, conc, npair = [], 0, [], [], [], [], 0, 0
        for n, pm in sp["per_molecule"].items():
            E = pm["e500" if stage == "500" else "econv"]
            rr = []
            for t in runs:
                if t == "p48":     # the stored-label trajectory: stored label (500) and this run's converged fit
                    rr.append(res[n]["label"] if stage == "500" else res[n]["labconv"])
                    continue
                z = np.load(CONV / f"traj_{t}" / "square_reg0.005" / f"{n}.npz")
                if stage == "conv" or int(z["nit"]) <= 500:
                    rr.append(float(z["resid"]))
                else:
                    rr.append(float(z["snap_resid"][list(z["snap_iters"]).index(500)]))
            rr, ee = np.array(rr), 100 * np.array([E[t] for t in runs])
            if np.ptp(rr) < 1e-9 or np.ptp(ee) < 1e-6:
                continue
            rhos.append(float(spearmanr(rr, ee)[0]))
            hit += int(np.argmin(rr) == np.argmax(ee))
            pick_lo.append(ee[np.argmin(rr)]); pick_hi.append(ee[np.argmax(rr)])
            mean.append(ee.mean()); best.append(ee.max())
            for i in range(len(rr)):
                for j in range(i + 1, len(rr)):
                    if abs(rr[i] - rr[j]) < 1e-6:
                        continue
                    npair += 1
                    conc += int((rr[i] < rr[j]) == (ee[i] > ee[j]))
        out["six_runs"][stage] = {
            "molecules": len(rhos), "spearman_resid_corr_mean": float(np.mean(rhos)),
            "spearman_resid_corr_median": float(np.median(rhos)), "spearman_negative": int(np.sum(np.array(rhos) < 0)),
            "lowest_resid_is_best_energy": hit, "pick_lowest_resid": float(np.mean(pick_lo)),
            "pick_highest_resid": float(np.mean(pick_hi)), "run_mean": float(np.mean(mean)),
            "oracle_best": float(np.mean(best)), "pairs": npair, "lower_resid_better_energy": conc / npair}

    # (b) label variants of the same molecule
    conc, dd = 0, []
    for n in ene:
        av = [k for k in LABS if k in ene[n] and k in res.get(n, {})]
        for i in range(len(av)):
            for j in range(i + 1, len(av)):
                dr, de = res[n][av[i]] - res[n][av[j]], ene[n][av[i]] - ene[n][av[j]]
                if abs(dr) < 1e-6 or abs(de) < 5e-7:        # identical labels (same parameters) are not a pair
                    continue
                conc += int((dr < 0) == (de > 0))
                dd.append((dr, de))
    dd = np.array(dd)
    out["label_variants"] = {"molecules": len(ene), "pairs": len(dd), "lower_resid_better_energy": conc / len(dd),
                             "spearman_dresid_dcorr": float(spearmanr(dd[:, 0], dd[:, 1])[0])}

    # (c) networks vs the converged label
    for c in NET:
        dR, dE = [], []
        for n in ene:
            if c in res.get(n, {}) and "labconv" in res[n] and c in ene[n]:
                dR.append(res[n][c] - res[n]["labconv"]); dE.append(100 * (ene[n][c] - ene[n]["labconv"]))
        dR, dE = np.array(dR), np.array(dE)
        out["networks_vs_labconv"][c] = {"molecules": len(dR), "resid_higher": int((dR > 0).sum()),
                                         "mean_dresid": float(dR.mean()), "energy_better": int((dE > 0).sum()),
                                         "mean_dcorr_points": float(dE.mean()),
                                         "resid_higher_and_energy_better": int(((dR > 0) & (dE > 0)).sum())}

    # (d) across molecules
    for s in sets:
        ns = [n for n in d[s]["per_molecule"] if "labconv" in d[s]["per_molecule"][n] and "labconv" in d[s]["resid"].get(n, {})]
        r = [d[s]["resid"][n]["labconv"] for n in ns]
        e = [d[s]["per_molecule"][n]["labconv"] for n in ns]
        out["across_molecules_labconv"][s] = {"n": len(ns), "spearman_resid_corr": float(spearmanr(r, e)[0])}

    print(json.dumps(out, indent=1))
    json.dump(out, open(NEW / "summary_resid_vs_energy.json", "w"), indent=1)
    print(f"-> {NEW / 'summary_resid_vs_energy.json'}")


if __name__ == "__main__":
    main()
