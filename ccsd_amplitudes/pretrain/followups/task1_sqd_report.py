#!/usr/bin/env python3
"""Tables for task 1 (noiseless QSCI with Lin et al.'s protocol), from the files written by task1_sqd_lin.py.

  python3 -m pretrain.followups.task1_sqd_report [--set smallval|n17|n17u] [--json out.json] [--md out.md]

Energies: LUCJ = exact variational energy (gpu_energy, complex64); QSCI = PySCF selected-CI energy in the sampled
subspace.  Protocols: "lin" (primary) 10^5 samples, 10 batches of 4,000, max_dim 4,000 (QSCI value of one draw =
mean over its 10 batches, her figures); "n2631g" 10^6 samples, one batch, no max_dim (her 16-orbital script);
"top100"/"top300": n2631g with max_dim = 100 / 300 strings per spin (fixed subspace size).  Seed 0 is her
entropy=0; seeds 1-4 are independent re-samplings ("5-seed" = mean over the draws that exist, sd = their spread:
5 draws for "lin"; for "n2631g" 3 draws (seeds 0-2) of the non-cdf_lin small-val states, else seed 0 only).
Set "n17u" = the norb-17 set without its three duplicate reactants (rxn4439/4440/4441/4442_R are one molecule).
Errors are E - E_ref in mHa.  "% corr" = 100 (E_HF - E) / (E_HF - E_ref) (share of the reference correlation
energy recovered); the LUCJ column uses the project's convention, % of the CCSD correlation energy.
"""
from __future__ import annotations

import argparse
import io
import json
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
TASK = "task1_sqd_vs_lin"
WORK = ROOT / "runs_ot" / "energy_tasks" / TASK
RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / TASK
CANDS = ["truncated", "cdf_lin", "label", "pre4", "rl4L", "rl4n29f"]
LAB = {"truncated": "truncated CCSD", "cdf_lin": "compressed DF (Lin)", "label": "optimize=True label",
       "pre4": "pretrained x4", "rl4L": "NN+RL rl4L", "rl4n29f": "NN+RL rl4n29f"}
SETS = {"smallval": "pretrain/rl/small_val.txt",
        "n17": "pretrain/opt_true/results/followups/task1_sqd_vs_lin/names_n17_c2hno.txt",
        "n17u": "pretrain/opt_true/results/followups/task1_sqd_vs_lin/names_n17_c2hno_unique.txt"}
NB = {"lin": 10}


def load_all():
    st = {}
    p = RES / "lucj_states.jsonl"
    if p.exists():
        for ln in open(p):
            r = json.loads(ln)
            st[(r["name"], r["cand"], r["seed"])] = r
    diag = {}
    for f in (WORK / "diag").glob("*.json"):
        name, cand, s, pr = f.stem.split("__")
        recs = json.load(open(f))["batches"]
        if len(recs) >= NB.get(pr, 1):                     # complete units only
            diag.setdefault((name, cand, int(s[1:])), {})[pr] = recs
    fci = json.load(open(RES / "fci_refs.json")) if (RES / "fci_refs.json").exists() else {}
    cct = json.load(open(ROOT / "pretrain" / "opt_true" / "results" / "ccsd_t_refs.json"))
    return st, diag, fci, cct


def ham_info(name):
    d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
    return dict(norb=int(d["norb"]), nelec=(int(d["nelec_a"]), int(d["nelec_b"])), e_hf=float(d["e_hf"]),
                e_ccsd=float(d["e_ccsd"]))


def sp(x, y):
    r = spearmanr(x, y)
    return float(getattr(r, "statistic", r[0]))


def build(names, st, diag, fci, cct):
    rows = {}
    for n in names:
        hi = ham_info(n)
        refs = dict(fci=fci.get(n, {}).get("e_fci"), ccsdt=cct.get(n, {}).get("e_ccsd_t"), ccsd=hi["e_ccsd"])
        for c in CANDS:
            seeds = sorted(s for (nn, cc, s) in st if nn == n and cc == c)
            if not seeds:
                continue
            s0 = st[(n, c, seeds[0])]
            r = dict(name=n, cand=c, norb=hi["norb"], nelec=hi["nelec"], e_hf=hi["e_hf"], refs=refs,
                     e_lucj=s0["e_lucj"], corr_lucj=100 * s0["corr_frac_ccsd"], p_hf=s0["p_hf"],
                     entropy=s0["entropy"], seeds={})
            for s in seeds:
                d = diag.get((n, c, s), {})
                rs = dict(n_unique_1e5=st[(n, c, s)]["protocols"]["lin"]["n_unique_subset"],
                          n_unique_1e6=st[(n, c, s)]["n_unique_1M"], t_energy=st[(n, c, s)]["t_energy"],
                          t_sample=st[(n, c, s)]["t_sample"], t_ci=st[(n, c, s)]["t_ci"])
                for pr, recs in d.items():
                    e = np.array([x["energy"] for x in recs])
                    rs[pr] = dict(e=e.tolist(), mean=float(e.mean()), min=float(e.min()), max=float(e.max()),
                                  dim=[x["dim_a"] for x in recs], n_distinct=len({x["sha1"] for x in recs}),
                                  s2max=float(max(x["spin_square"] for x in recs)),
                                  t=[x["t_solve"] for x in recs if not x.get("dedup")],
                                  threads=[x.get("threads") for x in recs if not x.get("dedup")])
                r["seeds"][s] = rs
            rows[(n, c)] = r
    return rows


def has(r, pr, seed=0):
    return r is not None and seed in r["seeds"] and pr in r["seeds"][seed]


def q_seeds(r, pr):
    """Per-seed QSCI values (batch means) of one (molecule, candidate)."""
    return np.array([r["seeds"][s][pr]["mean"] for s in sorted(r["seeds"]) if pr in r["seeds"][s]])


def err(e, ref):
    return None if (e is None or ref is None) else 1e3 * (e - ref)


def pct(e, e_hf, ref):
    return None if (e is None or ref is None) else 100 * (e_hf - e) / (e_hf - ref)


def fmt(x, f="{:.2f}"):
    return "–" if x is None or (isinstance(x, float) and np.isnan(x)) else f.format(x)


def table_per_molecule(rows, names, ref="fci", pr="lin"):
    print(f"| molecule | candidate | LUCJ % CCSD corr | LUCJ err | QSCI err seed 0: mean (min / max) | "
          f"distinct batches | 5-seed mean ± sd | QSCI % corr | dim/spin | p_HF |")
    print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for n in names:
        for c in CANDS:
            r = rows.get((n, c))
            if not has(r, pr):
                continue
            R = r["refs"][ref]
            q = r["seeds"][0][pr]
            ms = q_seeds(r, pr)
            sd = (f"{fmt(err(ms.mean(), R))} ± {1e3 * ms.std(ddof=1):.2f} ({len(ms)})" if len(ms) > 1 else "–")
            print(f"| {n} | {LAB[c]} | {r['corr_lucj']:.1f} | {fmt(err(r['e_lucj'], R), '{:.1f}')} | "
                  f"{fmt(err(q['mean'], R))} ({fmt(err(q['min'], R))} / {fmt(err(q['max'], R))}) | "
                  f"{q['n_distinct']} | {sd} | {fmt(pct(ms.mean(), r['e_hf'], R), '{:.2f}')} | "
                  f"{np.mean(q['dim']):.0f} | {r['p_hf']:.3f} |")


def matrix(rows, names, ref="fci", pr="lin"):
    """Compact per-molecule table: QSCI error (mean of the draws) with the LUCJ % CCSD corr in parentheses."""
    print("| molecule | " + " | ".join(LAB[c] for c in CANDS) + " | best QSCI | best LUCJ |")
    print("|:--|" + "--:|" * len(CANDS) + ":--|:--|")
    for n in names:
        cells, q, l = [], {}, {}
        for c in CANDS:
            r = rows.get((n, c))
            if not has(r, pr) or r["refs"][ref] is None:
                cells.append("–")
                continue
            q[c] = err(q_seeds(r, pr).mean(), r["refs"][ref])
            l[c] = r["corr_lucj"]
            cells.append(f"{q[c]:.1f} ({l[c]:.0f})")
        bq = LAB[min(q, key=q.get)] if q else "–"
        bl = LAB[max(l, key=l.get)] if l else "–"
        print(f"| {n} | " + " | ".join(cells) + f" | {bq} | {bl} |")


def summary(rows, names, ref="fci", pr="lin"):
    print(f"| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | "
          f"QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    lab = {n: q_seeds(rows[(n, "label")], pr).mean() for n in names if has(rows.get((n, "label")), pr)}
    for c in CANDS:
        rs = [rows[(n, c)] for n in names if has(rows.get((n, c)), pr) and rows[(n, c)]["refs"][ref] is not None]
        if not rs:
            continue
        el = [err(r["e_lucj"], r["refs"][ref]) for r in rs]
        e0 = [err(r["seeds"][0][pr]["mean"], r["refs"][ref]) for r in rs]
        e5 = [err(q_seeds(r, pr).mean(), r["refs"][ref]) for r in rs]
        sd = [1e3 * q_seeds(r, pr).std(ddof=1) for r in rs if len(q_seeds(r, pr)) > 1]
        pc = [pct(q_seeds(r, pr).mean(), r["e_hf"], r["refs"][ref]) for r in rs]
        dim = [np.mean(r["seeds"][0][pr]["dim"]) for r in rs]
        nu = [r["seeds"][0]["n_unique_1e5"] for r in rs]
        win = [q_seeds(r, pr).mean() < lab[r["name"]] for r in rs if r["name"] in lab]
        print(f"| {LAB[c]} | {len(rs)} | {np.mean([r['corr_lucj'] for r in rs]):.1f} | {np.mean(el):.1f} | "
              f"{np.mean(e0):.2f} | {np.mean(e5):.2f} | {fmt(np.mean(sd) if sd else None)} | {np.mean(pc):.2f} | "
              f"{np.mean(dim):.0f} | {np.mean([r['p_hf'] for r in rs]):.3f} | {np.mean(nu):.0f} | "
              f"{'—' if c == 'label' else f'{sum(win)}/{len(win)}'} |")


def pairwise(rows, names, pr="lin", base="label"):
    """Candidate vs base per molecule: dE_LUCJ and dE_QSCI (mHa, negative = candidate lower), with the seed noise."""
    print(f"| candidate vs {LAB[base]} | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | "
          f"QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|")
    for c in CANDS:
        if c == base:
            continue
        dl, dq, sig_lo, sig_hi = [], [], 0, 0
        for n in names:
            a, b = rows.get((n, base)), rows.get((n, c))
            if not (has(a, pr) and has(b, pr)):
                continue
            qa, qb = q_seeds(a, pr), q_seeds(b, pr)
            dl.append(1e3 * (b["e_lucj"] - a["e_lucj"]))
            d = 1e3 * (qb.mean() - qa.mean())
            dq.append(d)
            se = 1e3 * np.hypot(qa.std(ddof=1) / np.sqrt(len(qa)) if len(qa) > 1 else 0,
                                qb.std(ddof=1) / np.sqrt(len(qb)) if len(qb) > 1 else 0)
            sig_lo += d < -2 * se
            sig_hi += d > 2 * se
        if dl:
            print(f"| {LAB[c]} | {len(dl)} | {np.mean(dl):+.1f} | {sum(x < 0 for x in dl)}/{len(dl)} | "
                  f"{np.mean(dq):+.2f} | {sum(x < 0 for x in dq)}/{len(dq)} | {sig_lo}/{len(dl)} | {sig_hi}/{len(dl)} |")


def ranking(rows, names, pr="lin", cands=CANDS):
    """Spearman between the candidate order by LUCJ energy and by QSCI energy (5-seed means)."""
    per, xs, ys, conc, tot, conc_sig, tot_sig = [], [], [], 0, 0, 0, 0
    for n in names:
        cs = [c for c in cands if has(rows.get((n, c)), pr)]
        if len(cs) < 3:
            continue
        el = np.array([rows[(n, c)]["e_lucj"] for c in cs])
        qq = [q_seeds(rows[(n, c)], pr) for c in cs]
        qs = np.array([q.mean() for q in qq])
        se = np.array([q.std(ddof=1) / np.sqrt(len(q)) if len(q) > 1 else 0.0 for q in qq])
        per.append(sp(el, qs))
        xs += list(el - el.mean())
        ys += list(qs - qs.mean())
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                same = bool(np.sign(el[i] - el[j]) == np.sign(qs[i] - qs[j]))
                conc += same
                tot += 1
                if abs(qs[i] - qs[j]) > 2 * np.hypot(se[i], se[j]):
                    conc_sig += same
                    tot_sig += 1
    if not per:
        return None
    return dict(n_mol=len(per), n_cand=len(cands), rho_mean=float(np.mean(per)), rho_median=float(np.median(per)),
                rho_min=float(np.min(per)), rho_max=float(np.max(per)), rho_pooled=sp(xs, ys),
                concordant=f"{conc}/{tot}", concordant_resolved=f"{conc_sig}/{tot_sig}",
                per=[round(x, 2) for x in per])


def drivers(rows, names, pr="lin", cands=CANDS):
    """Pooled (within-molecule centred) Spearman of the QSCI energy with dim, p_HF, entropy and LUCJ energy."""
    q, d, p, h, e = [], [], [], [], []
    for n in names:
        cs = [c for c in cands if has(rows.get((n, c)), pr)]
        if len(cs) < 3:
            continue
        def cen(v):
            v = np.array(v, dtype=float)
            return list(v - v.mean())
        q += cen([q_seeds(rows[(n, c)], pr).mean() for c in cs])
        d += cen([np.mean([np.mean(rows[(n, c)]["seeds"][s][pr]["dim"]) for s in rows[(n, c)]["seeds"]
                           if pr in rows[(n, c)]["seeds"][s]]) for c in cs])
        p += cen([rows[(n, c)]["p_hf"] for c in cs])
        h += cen([rows[(n, c)]["entropy"] for c in cs])
        e += cen([rows[(n, c)]["e_lucj"] for c in cs])
    if not q:
        return None
    return dict(n=len(q), vs_dim=round(sp(q, d), 3), vs_p_hf=round(sp(q, p), 3), vs_entropy=round(sp(q, h), 3),
                vs_e_lucj=round(sp(q, e), 3))


def improvement(rows, names, ref="fci", pr="lin", base="truncated"):
    """Lin et al.'s improvement ratio: QSCI error of the baseline / QSCI error of the candidate (5-seed means)."""
    out = {}
    for c in CANDS:
        rat = []
        for n in names:
            a, b = rows.get((n, base)), rows.get((n, c))
            if not (has(a, pr) and has(b, pr)) or a["refs"][ref] is None:
                continue
            R = a["refs"][ref]
            ea, eb = q_seeds(a, pr).mean() - R, q_seeds(b, pr).mean() - R
            if ea > 0 and eb > 0:                        # ratios are meaningless for errors <= 0 (CCSD(T) ref)
                rat.append(ea / eb)
        if rat:
            out[c] = dict(n=len(rat), mean=float(np.mean(rat)), geomean=float(np.exp(np.mean(np.log(rat)))),
                          min=float(np.min(rat)), max=float(np.max(rat)))
    return out


def costs(rows, names):
    t = {"t_energy": [], "t_sample": [], "t_ci": []}
    td = {}
    for n in names:
        for c in CANDS:
            r = rows.get((n, c))
            if r is None:
                continue
            for s, rs in r["seeds"].items():
                for k in t:
                    t[k].append(rs[k])
                for pr in ("lin", "n2631g", "top100", "top300"):
                    if pr in rs:
                        for x, th in zip(rs[pr]["t"], rs[pr]["threads"]):
                            td.setdefault((pr, r["norb"], c == "cdf_lin"), []).append((x, th))
    print("| stage | median s | n |")
    print("|:--|--:|--:|")
    for k, v in t.items():
        if v:
            print(f"| {k} (per state / per draw) | {np.median(v):.1f} | {len(v)} |")
    for (pr, norb, cdf), v in sorted(td.items()):
        x = np.array([a for a, _ in v])
        th = sorted({b for _, b in v})
        print(f"| diag {pr}, norb {norb}{', cdf_lin' if cdf else ''} (threads {th}) | {np.median(x):.1f} "
              f"(max {x.max():.0f}) | {len(x)} |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="smallval", choices=list(SETS))
    ap.add_argument("--json", default=None)
    ap.add_argument("--md", default=None)
    a = ap.parse_args()
    names = [ln.strip() for ln in open(ROOT / SETS[a.set]) if ln.strip()]
    st, diag, fci, cct = load_all()
    rows = build(names, st, diag, fci, cct)
    refs = ["fci", "ccsdt"] if a.set == "smallval" else ["ccsdt"]
    res = {}
    buf = io.StringIO()
    with redirect_stdout(buf):
        for pr in ["lin", "n2631g", "top100", "top300"]:
            for ref in refs:
                print(f"\n### {a.set}: protocol {pr}, reference {ref}\n")
                summary(rows, names, ref, pr)
                imp = improvement(rows, names, ref, pr)
                if imp and pr in ("lin", "n2631g"):
                    print(f"\nQSCI error ratio truncated / candidate (5-seed means, reference {ref}):\n")
                    print("| candidate | n | mean | geometric mean | min | max |")
                    print("|:--|--:|--:|--:|--:|--:|")
                    for c, v in imp.items():
                        if c != "truncated":
                            print(f"| {LAB[c]} | {v['n']} | {v['mean']:.1f} | {v['geomean']:.1f} | {v['min']:.1f} | "
                                  f"{v['max']:.1f} |")
                res[f"{pr}_{ref}_improvement"] = imp
                if pr in ("lin", "n2631g"):
                    for base in ("label", "cdf_lin"):
                        ib = improvement(rows, names, ref, pr, base=base)
                        print(f"\nQSCI error ratio {LAB[base]} / candidate (reference {ref}; > 1 = candidate "
                              f"better; molecules with an error <= 0 skipped):\n")
                        print("| candidate | n | mean | geometric mean | min | max |")
                        print("|:--|--:|--:|--:|--:|--:|")
                        for c, v in ib.items():
                            if c not in (base, "truncated"):
                                print(f"| {LAB[c]} | {v['n']} | {v['mean']:.2f} | {v['geomean']:.2f} | "
                                      f"{v['min']:.2f} | {v['max']:.2f} |")
                        res[f"{pr}_{ref}_ratio_vs_{base}"] = ib
            print(f"\n### {a.set}: protocol {pr}, candidate vs label (QSCI: 5-seed means)\n")
            pairwise(rows, names, pr)
            print()
            for tag, cs in [("all 6", CANDS), ("without cdf_lin", [c for c in CANDS if c != "cdf_lin"]),
                            ("label + network family", ["label", "pre4", "rl4L", "rl4n29f"])]:
                rk = ranking(rows, names, pr, cs)
                dr = drivers(rows, names, pr, cs)
                if rk:
                    print(f"* ranking LUCJ vs QSCI ({tag}): rho per molecule mean {rk['rho_mean']:.2f} (median "
                          f"{rk['rho_median']:.2f}, min {rk['rho_min']:.2f}, max {rk['rho_max']:.2f}), pooled "
                          f"{rk['rho_pooled']:.2f}; concordant pairs {rk['concordant']} (resolved beyond 2 s.e.: "
                          f"{rk['concordant_resolved']}); per molecule {rk['per']}")
                if dr:
                    print(f"  pooled Spearman of QSCI energy with dim/spin {dr['vs_dim']}, p_HF {dr['vs_p_hf']}, "
                          f"entropy {dr['vs_entropy']}, LUCJ energy {dr['vs_e_lucj']} (n={dr['n']})")
                res[f"{pr}_ranking_{tag}"] = rk
                res[f"{pr}_drivers_{tag}"] = dr
        print(f"\n### {a.set}: per molecule, QSCI error in mHa (protocol lin, reference {refs[0]}, mean of the "
              f"draws) and LUCJ % CCSD corr in parentheses\n")
        matrix(rows, names, refs[0], "lin")
        print(f"\n### {a.set}: per molecule and candidate (protocol lin, reference {refs[0]}, mHa)\n")
        table_per_molecule(rows, names, refs[0], "lin")
        print(f"\n### {a.set}: costs\n")
        costs(rows, names)
    txt = buf.getvalue()
    print(txt)
    # every batch record of this set, copied out of the (gitignored) work dir into the results dir
    keep = ("batch", "energy", "dim_a", "subspace_dim", "spin_square", "c_hf", "t_solve", "t_s2", "threads", "host",
            "dedup", "sha1")
    allb = {f"{n}|{c}|s{s}|{pr}": [{k: x.get(k) for k in keep} for x in recs]
            for (n, c, s), d in diag.items() if n in names for pr, recs in d.items()}
    if a.set != "n17u":                                  # n17u is a subset of n17, whose batch file has every unit
        json.dump(dict(protocols=dict(lin="10^5 samples, 10 batches x 4000, max_dim 4000, 1 iteration",
                                      n2631g="10^6 samples, 1 batch, max_dim None, 1 iteration",
                                      top100="n2631g with max_dim 100", top300="n2631g with max_dim 300"),
                       solver="qiskit_addon_sqd.fermion.solve_sci(spin_sq=0.0), PySCF selected CI",
                       units=allb), open(RES / f"qsci_batches_{a.set}.json", "w"), indent=0)
    if a.md:
        Path(a.md).write_text(txt)
    if a.json:
        json.dump(dict(results=res, rows={f"{k[0]}|{k[1]}": v for k, v in rows.items()}), open(a.json, "w"),
                  indent=1, default=float)


if __name__ == "__main__":
    main()
