#!/usr/bin/env python3
"""Tables of docs/followups/task2_direct_opt.md from the task-2 result files (prints markdown, writes summary.json).

  python -m pretrain.followups.task2_report
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.followups import task2_common as C  # noqa: E402

TRACE_AT = (1, 50, 100, 200, 300, 500)


def load(obj):
    out = {}
    for f in sorted((C.RESULTS / obj).glob("*.json")):
        r = json.load(open(f))
        tag = r["args"].get("tag") or ""
        out[(r["name"], r["start"], r["optimizer"], tag)] = r
    return out


_FCI = None


def fci_e(name):
    """FCI energy from Task 1 (pyscf direct_spin0, all nine norb 15-16 molecules), or None."""
    global _FCI
    if _FCI is None:
        p = C.ROOT / "pretrain/opt_true/results/followups/task1_sqd_vs_lin/fci_refs.json"
        _FCI = json.load(open(p)) if p.exists() else {}
    v = _FCI.get(name)
    return None if v is None else v["e_fci"]


def mha(E, name):
    ef = fci_e(name)
    return None if ef is None else 1000 * (E - ef)


def pct(E, r):
    return C.corr_pct(E, r["e_hf"], r["e_ccsd"])


def fmt(x, nd=2, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def table_lucj(res):
    """Objective A.  SPSA columns: 'cal' = gain calibrated at the start point (first step 0.008 per coordinate),
    'aL' = gain calibrated at the label start of the same molecule (tag aLabel; identical to 'cal' for the label)."""
    names = sorted({k[0] for k in res if not k[3]}, key=lambda n: (C.load_ham(n)["norb"], n))
    print("\n### Objective A: exact LUCJ energy (% of CCSD correlation energy), 500 evaluations\n")
    print("| molecule | norb | start | start % | NOMAD | Δ | SPSA (cal.) | Δ | SPSA (label gain) | Δ |")
    print("|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|")
    rows = []

    def get(n, st, o, tag=""):
        return res.get((n, st, o, tag))

    for n in names:
        for st in ("label", "rl4n29f"):
            a, b = get(n, st, "nomad"), get(n, st, "spsa")
            c = b if st == "label" else get(n, st, "spsa", "aLabel")
            if not (a or b or c):
                continue
            r0 = a or b or c
            c0 = r0["corr_var0"]
            v = [x["corr_var_best"] if x else None for x in (a, b, c)]
            rows.append(dict(name=n, start=st, c0=c0, nomad=v[0], spsa=v[1], spsa_aL=v[2]))
            cells = " | ".join(f"{fmt(x)} | {fmt(None if x is None else x - c0, sign=True)}" for x in v)
            print(f"| {n} | {r0['norb']} | {st} | {c0:.2f} | {cells} |")
    full = [r["name"] for r in rows if r["start"] == "rl4n29f" and None not in (r["nomad"], r["spsa_aL"])
            and any(q["name"] == r["name"] and q["start"] == "label" and None not in (q["nomad"], q["spsa"])
                    for q in rows)]
    if full:
        print(f"\nMeans over the {len(full)} molecules with NOMAD and SPSA from both starts:\n")
        print("| start | one call / start | NOMAD | SPSA (cal.) | SPSA (label gain) | best of the runs |")
        print("|:--|--:|--:|--:|--:|--:|")
        for st in ("label", "rl4n29f"):
            rr = [r for r in rows if r["start"] == st and r["name"] in full]
            c0 = np.mean([r["c0"] for r in rr])
            cols = []
            for key in ("nomad", "spsa", "spsa_aL"):
                vals = [r[key] for r in rr]
                cols.append("–" if None in vals else f"{np.mean(vals):.2f} ({np.mean(vals) - c0:+.2f})")
            best = np.mean([max(x for x in (r["nomad"], r["spsa"], r["spsa_aL"]) if x is not None) for r in rr])
            print(f"| {st} | {c0:.2f} | " + " | ".join(cols) + f" | {best:.2f} ({best - c0:+.2f}) |")
    # best-so-far traces
    print("\nBest-so-far (mean % corr over molecules) after k evaluations:\n")
    print("| start / optimizer | " + " | ".join(f"k={k}" for k in TRACE_AT) + " |")
    print("|:--|" + "--:|" * len(TRACE_AT))
    for st, o, tag in (("label", "nomad", ""), ("label", "spsa", ""), ("rl4n29f", "nomad", ""),
                       ("rl4n29f", "spsa", ""), ("rl4n29f", "spsa", "aLabel")):
        rs = [res[(n, st, o, tag)] for n in names if (n, st, o, tag) in res]
        if not rs:
            continue
        vals = [np.mean([pct(r["best_trace"][min(k, len(r["best_trace"])) - 1], r) for r in rs]) for k in TRACE_AT]
        lab = f"{st} / {o}" + (f" ({tag})" if tag else "")
        print(f"| {lab} ({len(rs)} mol.) | " + " | ".join(f"{v:.2f}" for v in vals) + " |")
    # long runs
    longs = {k: r for k, r in res.items() if k[3].startswith("long")}
    if longs:
        ks = (1, 500, 1000, 2000, 3000, 5000)
        print("\nLong SPSA runs (best-so-far % corr after k evaluations):\n")
        print("| run | " + " | ".join(f"k={k}" for k in ks) + " |")
        print("|:--|" + "--:|" * len(ks))
        for k_, r in sorted(longs.items()):
            tr = r["best_trace"]
            vals = [fmt(pct(tr[min(k, len(tr)) - 1], r)) if len(tr) >= min(k, len(tr)) else "–" for k in ks]
            print(f"| {k_[0]} {k_[1]} {k_[3]} ({r['n_evals']} evals, {r['wall_s'] / 3600:.1f} h) | " + " | ".join(vals) + " |")
    return rows


def table_match(res):
    """Objective A: how many evaluations a per-molecule optimizer from the label needs to reach the one-call value."""
    names = sorted({k[0] for k in res if k[1] == "label" and not k[3]}, key=lambda n: (C.load_ham(n)["norb"], n))
    print("\n### Evaluations needed from the label to match one network call (objective A, 500-evaluation runs)\n")
    print("| molecule | label % | one call % | NOMAD best % | NOMAD: evals to match | SPSA best % | "
          "SPSA: evals to match (min) |")
    print("|:--|--:|--:|--:|--:|--:|--:|")
    hits = {"nomad": [], "spsa": []}
    for n in names:
        one = next((res[(n, "rl4n29f", o, t)]["corr_var0"] for o, t in (("nomad", ""), ("spsa", ""), ("spsa", "aLabel"))
                    if (n, "rl4n29f", o, t) in res), None)
        if one is None:
            continue
        cells, lab0 = [], None
        for o in ("nomad", "spsa"):
            r = res.get((n, "label", o, ""))
            if r is None:
                cells += ["–", "–"]
                continue
            lab0 = r["corr_var0"]
            tr = np.array(pct_trace_r(r))
            idx = np.nonzero(tr >= one)[0]
            if len(idx):
                k = int(idx[0]) + 1
                hits[o].append((k, k * r["wall_s"] / r["n_evals"] / 60))
                cell = f"{k}" if o == "nomad" else f"{k} ({k * r['wall_s'] / r['n_evals'] / 60:.1f})"
            else:
                hits[o].append(None)
                cell = f"not in {r['n_evals']}"
            cells += [f"{tr.max():.2f}", cell]
        print(f"| {n} | {lab0:.2f} | {one:.2f} | " + " | ".join(cells) + " |")
    for o, h in hits.items():
        ok = [x for x in h if x is not None]
        if h:
            med = f"; median {np.median([x[0] for x in ok]):.0f} evaluations ({np.median([x[1] for x in ok]):.1f} min)" if ok else ""
            print(f"\n* {o.upper() if o == 'nomad' else 'SPSA'} from the label reached the one-call value on "
                  f"{len(ok)}/{len(h)} molecules{med}.")


def pct_trace_r(r):
    return [C.corr_pct(v, r["e_hf"], r["e_ccsd"]) for v in r["best_trace"]]


def table_tune(res):
    tun = {k: r for k, r in res.items() if k[3].startswith("tune")}
    if not tun:
        return
    print("\n### SPSA tuning (training molecule, label start, objective A, 500 evaluations)\n")
    print("| molecule | c | first-step target | start | best | Δ |")
    print("|:--|--:|--:|--:|--:|--:|")
    for k, r in sorted(tun.items(), key=lambda kv: (kv[1]["spsa"]["c"], kv[1]["spsa"]["target"])):
        s = r["spsa"]
        print(f"| {k[0]} | {s['c']} | {s['target']} | {r['corr_var0']:.2f} | {r['corr_var_best']:.2f} | "
              f"{r['corr_var_best'] - r['corr_var0']:+.2f} |")


def table_qsci(res):
    if not res:
        return
    print("\n### Objective B: QSCI energy of 10^4 exact samples (her TN-optimization objective)\n")
    print("| molecule | start | optimizer | evals | wall (h) | s / eval | QSCI start → best (% corr) | "
          "QSCI − FCI start → best (mHa) | best at | variational % start → at best | strings per spin start → at best "
          "(max) |")
    print("|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for k, r in sorted(res.items()):
        if k[3]:
            continue
        hist = [json.loads(ln) for ln in open(C.LOGS / "qsci" / f"{k[0]}__{k[1]}__{k[2]}.jsonl")]
        dmax = max(h["dim"][0] for h in hist)
        note = " (stopped)" if (r.get("stopped_early") or r.get("truncated")) else ""
        print(f"| {k[0]} | {k[1]} | {k[2]} | {r['n_evals']}{note} | {r['wall_s'] / 3600:.2f} | {r['t_eval_mean']:.1f} | "
              f"{r['corr_f0']:.2f} → {r['corr_f_best']:.2f} | {fmt(mha(r['f0'], k[0]))} → {fmt(mha(r['f_best'], k[0]))} | "
              f"{r['best_at']} | {r['corr_var0']:.2f} → {r['corr_var_best']:.2f} | {hist[0]['dim'][0]} → "
              f"{hist[r['best_at'] - 1]['dim'][0]} ({dmax}) |")


def table_scores():
    p = C.RESULTS / "score_task1_protocol.jsonl"
    if not p.exists():
        return {}
    rows = {}
    for ln in open(p):
        r = json.loads(ln)
        rows[(r["name"], r["cand"])] = r
    names = sorted({k[0] for k in rows}, key=lambda n: (C.load_ham(n)["norb"], n))
    order = ["start:label", "lucj:label:nomad", "lucj:label:spsa", "qsci:label:nomad", "qsci:label:spsa",
             "start:rl4n29f", "lucj:rl4n29f:nomad", "lucj:rl4n29f:spsa", "qsci:rl4n29f:nomad", "qsci:rl4n29f:spsa"]
    print("\n### Task-1 protocol (10^5 exact samples, 10 batches x 4000, her SQD settings)\n")
    print("| molecule | parameters | variational % | SQD mean % [range over batches] | SQD mean − FCI (mHa) | "
          "SQD mean − CCSD(T) (mHa) | unique / 10^5 | strings |")
    print("|:--|:--|--:|--:|--:|--:|--:|--:|")
    for n in names:
        for c in order + sorted({k[1] for k in rows if k[0] == n} - set(order)):
            r = rows.get((n, c))
            if not r:
                continue
            lo, hi = sorted((r['corr_sqd_min'], r['corr_sqd_max']))
            print(f"| {n} | {c} | {r['corr_var']:.2f} | {r['corr_sqd_mean']:.2f} [{lo:.2f}, "
                  f"{hi:.2f}] | {fmt(mha(r['E_sqd_mean'], n))} | {fmt(r['err_sqd_mean_vs_ccsdt_mHa'], 1)} | "
                  f"{r['n_unique_100k']} | {min(r['dims'])}–{max(r['dims'])} |")
    # means over the molecules that have every objective-A parameter set
    sets = ["start:label", "lucj:label:spsa", "start:rl4n29f", "lucj:rl4n29f:spsa:aLabel"]
    full = [n for n in names if all((n, c) in rows for c in sets)]
    if full:
        print(f"\nMeans over the {len(full)} molecules scored for all four of {', '.join(sets)}:\n")
        print("| parameters | variational % | SQD mean − FCI (mHa) | strings per spin |")
        print("|:--|--:|--:|--:|")
        for c in sets + ["lucj:label:nomad", "lucj:rl4n29f:nomad"]:
            ns = [n for n in full if (n, c) in rows]
            if not ns:
                continue
            rr = [rows[(n, c)] for n in ns]
            print(f"| {c} ({len(ns)} mol.) | {np.mean([r['corr_var'] for r in rr]):.2f} | "
                  f"{np.mean([mha(r['E_sqd_mean'], r['name']) for r in rr]):.2f} | "
                  f"{np.mean([np.mean(r['dims']) for r in rr]):.0f} |")
    return rows


def table_dimscan():
    p = C.RESULTS / "dimscan.jsonl"
    if not p.exists():
        return
    rows = {}
    for ln in open(p):
        r = json.loads(ln)
        rows[(r["name"], r["cand"], r["k"])] = r
    names = sorted({k[0] for k in rows}, key=lambda n: (C.load_ham(n)["norb"], n))
    print("\n### Subspace energy at fixed size: k most probable strings (exact marginals), k x k subspace\n")
    for n in names:
        ks = sorted({k[2] for k in rows if k[0] == n})
        cands = sorted({k[1] for k in rows if k[0] == n}, key=lambda c: (not c.startswith("start"), c))
        has_fci = any(rows[k]["err_fci_mHa"] is not None for k in rows if k[0] == n)
        unit = "error vs FCI (mHa)" if has_fci else "% CCSD corr."
        print(f"\n{n} ({unit}; variational % in the second column):\n")
        print("| parameters | variational % | " + " | ".join(f"k={k}" for k in ks) + " |")
        print("|:--|--:|" + "--:|" * len(ks))
        for c in cands:
            vals = []
            for k in ks:
                r = rows.get((n, c, k))
                if not r:
                    vals.append("–")
                else:
                    vals.append(f"{r['err_fci_mHa']:.2f}" if has_fci else f"{r['corr']:.2f}")
            rv = next(rows[(n, c, k)] for k in ks if (n, c, k) in rows)
            print(f"| {c} | {rv['corr_var']:.2f} | " + " | ".join(vals) + " |")


def table_amortization(lucj, qsci):
    print("\n### Amortization: one network call vs per-molecule optimization\n")
    inf = {}
    for d in ("gpu", "cpu"):
        f = C.RESULTS / f"infer_time_{d}.json"
        if f.exists():
            inf[d] = json.load(open(f))
    for d, v in inf.items():
        print(f"* rl4n29f inference ({v['device']}): median {1000 * v['median_total_s']:.0f} ms per molecule "
              f"(input build incl. 3 frozen recycles + final recycle + realize); model load {v['t_load_model_s']:.1f} s.")
    for obj, res in (("A (LUCJ)", lucj), ("B (QSCI)", qsci)):
        for o in ("nomad", "spsa"):
            for st in ("label", "rl4n29f"):
                rs = [r for k, r in res.items() if k[1] == st and k[2] == o and k[3] in ("", "aLabel")]
                if not rs:
                    continue
                for norb in sorted({r["norb"] for r in rs}):
                    ws = [r["wall_s"] for r in rs if r["norb"] == norb]
                    te = [r["t_eval_mean"] for r in rs if r["norb"] == norb]
                    print(f"* objective {obj}, {o}, start {st}, norb {norb} ({len(ws)} run{'s' if len(ws) != 1 else ''}): "
                          f"{np.mean(ws) / 60:.1f} min per run ({np.mean(te):.2f} s per evaluation)")


def main():
    lucj = load("lucj")
    qsci = load("qsci")
    rows_a = table_lucj(lucj)
    table_match(lucj)
    table_tune(lucj)
    table_qsci(qsci)
    scores = table_scores()
    table_dimscan()
    table_amortization(lucj, qsci)
    summ = dict(lucj={"|".join(k): {kk: v for kk, v in r.items() if kk != "best_trace"} for k, r in lucj.items()},
                qsci={"|".join(k): {kk: v for kk, v in r.items() if kk != "best_trace"} for k, r in qsci.items()},
                lucj_rows=rows_a, scores=["|".join(k) for k in scores])
    C.dumpj(summ, C.RESULTS / "summary.json")


if __name__ == "__main__":
    main()
