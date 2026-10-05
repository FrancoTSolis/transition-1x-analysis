#!/usr/bin/env python3
"""Tables of docs/followups/sampler.md from pretrain/opt_true/results/followups/sampler/*.json.

  python3 -m pretrain.followups.sampler_report
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "sampler"


def short(name):
    return name.replace("_rxn", " ")


def table_validation(recs):
    print("| molecule (norb) | χ | E − E_exact (mHa) | discarded_sum | 1 − F | TVD full | TVD top-10³ (exact p) | "
          "TVD top-10³ (samples) | unique configs | unique α∪β | state build (s) | sampling 10⁵ (s) |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for r in recs:
        ex = r["exact_chi"]
        floor = ex["tvd_topK_sample_floor"]["1000"]
        u1 = ex["unique"]["1"]
        print(f"| {short(r['name'])} ({r['norb']}) | exact | 0 | – | – | – | – | {floor:.4f} (2nd exact set) | "
              f"{u1['unique']} / {ex['unique']['2']['unique']} | {u1['unique_union']} / {ex['unique']['2']['unique_union']} | – | – |")
        for chi, m in r["mps"].items():
            onemf = 1 - m["fidelity_from_exact_samples"]
            print(f"| | {chi} | {m['E_minus_exact_mHa']:+.3f} | {m['discarded_sum']:.1e} | {onemf:.1e} ± "
                  f"{2 * m['overlap_abs_se']:.0e} | {m['tvd_full_mc']:.1e} | {m['tvd_topK_exact']['1000']:.1e} | "
                  f"{m['tvd_topK_samples_vs_exact_p']['1000']:.4f} | {m['unique']['unique']} | {m['unique']['unique_union']} | "
                  f"{m['t_state']:.0f} | {m['t_sample']:.2f} |")


def table_sqd(recs, sqds):
    sets = ["exact1", "exact2", "mps64", "mps128", "mps256", "mps512", "mo"]
    print("| molecule | E_LUCJ − E_CCSD(T) | " + " | ".join(
        {"exact1": "exact χ #1", "exact2": "exact χ #2", "mo": "exact MO"}.get(s, f"MPS χ{s[3:]}") for s in sets)
        + " |")
    print("|:--|--:|" + "--:|" * len(sets))
    for r in recs:
        q = sqds.get(r["name"])
        if q is None:
            continue
        ref = r.get("e_ccsd_t") or r["e_ccsd"]
        row = [f"{1e3 * (r['exact_mo']['E'] - ref):+.2f}"]
        for s in sets:
            x = q["sets"].get(s)
            if x is None:
                row.append("–")
                continue
            row.append(f"{1e3 * (x['E_min'] - ref):+.3f} ({x['batches'][0]['dim']:,})")
        print(f"| {short(r['name'])} ({r['norb']}) | " + " | ".join(row) + " |")


def pooled_sets(name, cand="rl4n29f"):
    """QSCI energies of every independent sample set per source: the seeds file (6 sets) + the first draw(s) of
    sqd_<name>__<cand>.json (exact1/exact2 -> exact, mps<chi> -> mps<chi>).  Returns {source: [(E, dim), ...]}."""
    out = {}
    files = [RES / f"sqd_{name}__{cand}.json"] + sorted(RES.glob(f"sqd_{name}__{cand}_seeds*.json"))
    for f in files:
        if not f.exists():
            continue
        for st, x in json.load(open(f))["sets"].items():
            src = st.split("_s")[0]
            src = "exact" if src in ("exact1", "exact2") else src
            if not x.get("batches"):
                continue
            out.setdefault(src, []).append((x["E_min"], x["batches"][0]["dim"]))
    return out


def table_seeds(recs):
    """Pooled QSCI energies per source (mean ± s.e.m. over independent sample sets), then the differences of each MPS
    χ and of the MO basis from exact χ-basis sampling (independent seeds, so s.e. = hypot of the two)."""
    srcs = ["exact", "mps64", "mps128", "mps256", "mps512", "mo"]
    name = {"exact": "exact χ", "mo": "exact MO"}
    print("| molecule | " + " | ".join(name.get(s, f"MPS χ{s[3:]}") for s in srcs) + " |")
    print("|:--|" + "--:|" * len(srcs))
    diffs = {s: [] for s in srcs[1:]}
    drows = []
    for r in recs:
        P = pooled_sets(r["name"], r["cand"])
        ref = r.get("e_ccsd_t") or r["e_ccsd"]
        row, stats = [], {}
        for s in srcs:
            v = P.get(s, [])
            if not v:
                row.append("–")
                continue
            E = 1e3 * (np.array([e for e, _ in v]) - ref)
            m = E.mean()
            se = E.std(ddof=1) / np.sqrt(len(E)) if len(E) > 1 else float("nan")
            stats[s] = (m, se, len(E))
            dim = np.mean([d for _, d in v])
            row.append(f"{m:+.2f} ± {se:.2f} (n={len(E)}, {dim / 1e3:.0f}k)" if len(E) > 1 else f"{m:+.2f} (n=1, {dim / 1e3:.0f}k)")
        print(f"| {short(r['name'])} ({r['norb']}) | " + " | ".join(row) + " |")
        drow = []
        for s in srcs[1:]:
            if s in stats and "exact" in stats:
                d = stats[s][0] - stats["exact"][0]
                e = float(np.hypot(stats[s][1], stats["exact"][1]))
                diffs[s].append((d, e))
                drow.append(f"{d:+.2f} ± {e:.2f}")
            else:
                drow.append("–")
        drows.append(f"| {short(r['name'])} ({r['norb']}) | " + " | ".join(drow) + " |")
    print("\nDifference from exact χ-basis sampling (mHa; mean ± s.e.; last row: mean over the molecules):\n")
    print("| molecule | " + " | ".join(f"{name.get(s, f'MPS χ{s[3:]}')} − exact χ" for s in srcs[1:]) + " |")
    print("|:--|" + "--:|" * (len(srcs) - 1))
    for x in drows:
        print(x)
    mean = []
    for s in srcs[1:]:
        v = diffs[s]
        if v:
            mean.append(f"{np.mean([d for d, _ in v]):+.2f} ± {np.sqrt(sum(e * e for _, e in v)) / len(v):.2f}")
        else:
            mean.append("–")
    print("| mean | " + " | ".join(mean) + " |")


def main():
    recs = [json.load(open(f)) for f in sorted(RES.glob("validate_*.json"))]
    sqds = {}
    for f in sorted(RES.glob("sqd_*.json")):
        if "maxdim" in f.stem or "_seeds" in f.stem:
            continue
        q = json.load(open(f))
        sqds[q["name"]] = q
    print("## validation\n")
    table_validation(recs)
    print("\n## QSCI (mHa vs CCSD(T); subspace dimension)\n")
    table_sqd(recs, sqds)
    print("\n## QSCI pooled over independent sample sets (mean ± s.e.m., mHa vs CCSD(T))\n")
    table_seeds(recs)
    print("\n## n29\n")
    table_n29()
    return 0


def table_n29():
    refs = {}
    f = RES / "ccsd_t_refs_sampler.json"
    if f.exists():
        refs = json.load(open(f))
    print("| molecule | χ | E (Ha) | E − E_CCSD(T) (mHa) | % CCSD corr | discarded_sum | state (s) | 10⁵ / 10⁶ samples (s) | "
          "1 − F vs χ512 | TVD vs χ512 | distinct configs | α∪β strings | QSCI subspace (each of 10 batches) |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for f in sorted(RES.glob("n29_*.json")):
        r = json.load(open(f))
        et = refs.get(r["name"], {}).get("e_ccsd_t")
        for chi in sorted(r["mps"], key=int):
            m = r["mps"][chi]
            dE = f"{1e3 * (m['E'] - et):+.1f}" if (et and "E" in m) else "–"
            fid = f"{1 - m['fidelity_vs_ref']:.1e} ± {2 * m['overlap_se']:.0e}" if m.get("fidelity_vs_ref") else "(ref)"
            tvd = f"{m['tvd_vs_ref']:.1e}" if m.get("tvd_vs_ref") else "(ref)"
            u = m["unique"]
            dims = sorted(set(m["sqd_dims"]))
            print(f"| {short(r['name'])} | {chi} | {m.get('E', float('nan')):.6f} | {dE} | {100 * m.get('corr_frac', float('nan')):.2f} | "
                  f"{m['discarded_sum']:.1e} | {m['t_state']:.0f} | {m['t_sample']:.2f} / {m.get('t_sample_throughput', float('nan')):.1f} | "
                  f"{fid} | {tvd} | {u['unique']:,} | {u['unique_union']:,} | {', '.join(f'{d:,}' for d in dims)} |")
    for f in sorted(RES.glob("sqd_*maxdim*.json")):
        q = json.load(open(f))
        for st, x in q["sets"].items():
            b = x["batches"][0]
            et = refs.get(q["name"], {}).get("e_ccsd_t")
            print(f"\nQSCI {q['name']} {st} max_dim {q['settings']['max_dim']}: E {b['E']:.6f}"
                  + (f" ({1e3 * (b['E'] - et):+.1f} mHa vs CCSD(T))" if et else "")
                  + f", dim {b['dim']:,}, {b['t']:.0f} s, S^2 {b['spin_sq']:.3f}")
    for f in sorted(RES.glob("sci_bench_*.json")):
        q = json.load(open(f))
        print(f"\nsigma bench {q['name']} {q['set']} ({q['threads']} threads), full batch strings {q['full_strings']}:")
        for row in q["rows"]:
            print(f"  max_dim {row['max_dim']}: dim {row['dim']:,}  t_sigma {row['t_sigma']:.1f} s  "
                  f"(N-2 intermediates {row['dd_inter_a']:,})")


if __name__ == "__main__":
    sys.exit(main())
