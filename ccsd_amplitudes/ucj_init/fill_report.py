"""docs/ucj_init/report_template.md -> docs/ucj_init/report.md (generated tables), and the summary statistics that the
report's answer section quotes (printed, to check the text against the data)."""
import io
import json
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

from ucj_init import report_small
from ucj_init.report_exact import ORDER, evals_to, final_table, load, paired_table

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "ucj_init"
RES = ROOT / "pretrain" / "opt_true" / "results" / "ucj_init"


def capacity_table():
    rows = json.load(open(RES / "ansatz_capacity.json"))
    c2 = RES / "ansatz_capacity_c2.json"
    if c2.exists():
        rows += json.load(open(c2))
    kinds = [("real1", "his ansatz: real U, 1 layer, spin-free"), ("cplx1", "complex U, 1 layer"),
             ("lucj2sq", "LUCJ, 2 layers, square (ours)"), ("ucj2", "UCJ, 2 layers, dense J")]
    by = {}
    for r in rows:
        by.setdefault(r["name"], {})[r["kind"]] = r["best_pct"]
    out = ["| system | " + " | ".join(k[1] for k in kinds) + " |", "|:--|" + "--:|" * len(kinds)]
    for n, d in by.items():
        out.append(f"| {n} | " + " | ".join(f"{d[k]:.1f}" if k in d else "(not run)" for k, _ in kinds) + " |")
    return "\n".join(out)


def avijit_table():
    buf = io.StringIO()
    with redirect_stdout(buf):
        report_small.main()
    txt = buf.getvalue().split("\n\n")[0]
    return txt.strip()


def stats():
    runs = load("lbfgs")
    names = sorted({n for n, _ in runs})
    ref = "rl4n29f"
    print("molecules with rl4n29f:", len([n for n in names if (n, ref) in runs]), names)
    print("rl4n29f start range: %.1f-%.1f" % (min(runs[(n, ref)]["best"][0] for n in names if (n, ref) in runs),
                                             max(runs[(n, ref)]["best"][0] for n in names if (n, ref) in runs)))
    for init in ORDER:
        mols = [n for n in names if (n, init) in runs and (n, ref) in runs]
        if not mols or init == ref:
            continue
        fin = [(runs[(n, ref)]["best"][-1], runs[(n, init)]["best"][-1]) for n in mols]
        wins = sum(a > b + 0.05 for a, b in fin)
        ties = sum(abs(a - b) <= 0.05 for a, b in fin)
        b100 = [(runs[(n, ref)]["best"][100], runs[(n, init)]["best"][100]) for n in mols
                if len(runs[(n, ref)]["best"]) > 100 and len(runs[(n, init)]["best"]) > 100]
        e80 = [evals_to(runs[(n, init)]["best"], 80) for n in mols]
        print(f"{init:10s} n={len(mols):2d}  final: net ahead {wins}, tie {ties}, behind {len(fin) - wins - ties} | "
              f"at 100: net ahead {sum(a > b for a, b in b100)}/{len(b100)} | evals to 80%: "
              f"{sorted(e if e is not None else 999 for e in e80)}")
    t = load("trn")
    for init in ORDER:
        mols = [n for n, i in t if i == init]
        if mols:
            print("trn", init, "evals to 70%:", [evals_to(t[(n, init)]["best"], 70) for n in mols],
                  "final:", [round(t[(n, init)]["best"][-1], 1) for n in mols])


def main():
    s = (DOC / "report_template.md").read_text()
    ts = lambda n: n.endswith("_TS")
    lb = ("**All molecules** (norb 15–18):\n\n" + paired_table("lbfgs", budgets=(0, 10, 25, 50, 100, 150)) +
          "\n\n**TS molecules, norb 15–16:**\n\n" +
          paired_table("lbfgs", budgets=(0, 10, 25, 50, 100, 200),
                       names_filter=lambda n: ts(n) and not n.startswith(("C2HNO", "C2H2O2"))) +
          "\n\n**TS molecules, norb 17–18:**\n\n" +
          paired_table("lbfgs", budgets=(0, 10, 25, 50, 100, 150),
                       names_filter=lambda n: ts(n) and n.startswith(("C2HNO", "C2H2O2"))))
    s = s.replace("FILL_LBFGS", lb)
    s = s.replace("FILL_TRN", "Norb 15 (C2H3N P, TS) and norb 16 (C3H4 TS, C2H4O_rxn0724 TS), 400 evaluations:\n\n" +
                  paired_table("trn", budgets=(0, 10, 25, 50, 100, 200, 300)))
    s = s.replace("FILL_JONLY", final_table("lbfgsJ"))
    s = s.replace("FILL_FINAL", final_table("lbfgs", inits=("zero", "random0", "truncated", "cdf_lin", "label",
                                                             "pre4", "rl4L", "rl4n29f")))
    s = s.replace("FILL_CAPACITY", capacity_table())
    s = s.replace("FILL_AVIJIT", avijit_table())
    (DOC / "report.md").write_text(s)
    print("wrote", DOC / "report.md")
    stats()


if __name__ == "__main__":
    main()
