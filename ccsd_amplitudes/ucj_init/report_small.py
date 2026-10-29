"""Markdown tables of the small-system studies (avijit_exact.json, ansatz_capacity.json, rq2_test.json)."""
import json
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / "pretrain" / "opt_true" / "results" / "ucj_init"


def main():
    f = RES / "avijit_exact.json"
    if f.exists():
        rows = json.load(open(f))
        seeds = []
        for r in rows:
            if r["seed"] not in seeds:
                seeds.append(r["seed"])
        print("| system | FCI corr. (Eh) | " + " | ".join(seeds) + " |")
        print("|:--|--:|" + "--:|" * len(seeds))
        by = {}
        for r in rows:
            by.setdefault(r["name"], {})[r["seed"]] = r
        for name, d in by.items():
            fci = next(iter(d.values()))["e_fci"] - next(iter(d.values()))["e_hf"]
            cells = []
            for s in seeds:
                r = d.get(s)
                if r is None:
                    cells.append("")
                    continue
                ex = r["exact"]
                pct = 100 * ex["e_corr"] / fci
                stat = "" if ex["grad_inf"] < 1e-6 else f" (grad {ex['grad_inf']:.0e})"
                cells.append(f"{pct:.0f} % / {ex['cycles']}{stat}")
            print(f"| {name} | {fci:+.3f} | " + " | ".join(cells) + " |")
        n_bch2_div = sum(1 for r in rows if not r["bch2"]["converged"])
        print(f"\nBCH2 (same starts, same driver): {n_bch2_div}/{len(rows)} runs not converged; "
              f"min E_BCH2 - E_HF over those: {min(r['bch2']['e_corr'] for r in rows):.1f} Eh; "
              f"exact energy at the BCH2 end points above E_HF in "
              f"{sum(1 for r in rows if r['bch2']['E_exact_end'] - r['e_hf'] > 0)}/{len(rows)} runs\n")
    f = RES / "ansatz_capacity.json"
    if f.exists():
        rows = json.load(open(f))
        kinds = ["real1", "cplx1", "lucj2sq", "ucj2"]
        by = {}
        for r in rows:
            by.setdefault(r["name"], {})[r["kind"]] = r["best_pct"]
        print("| system | " + " | ".join(kinds) + " |\n|:--|" + "--:|" * len(kinds))
        for name, d in by.items():
            print(f"| {name} | " + " | ".join(f"{d[k]:.1f}" if k in d else "" for k in kinds) + " |")
    f = RES / "rq2_test.json"
    if f.exists():
        rows = json.load(open(f))
        print("\n| system | CISD | RQ2 minima | exact energy at those parameters |\n|:--|--:|:--|:--|")
        for r in rows:
            print(f"| {r['name']} | {r['cisd_pct']:.1f} | " + ", ".join(f"{x['pct']['rq2']:.0f}" for x in r["runs"]) +
                  " | " + ", ".join(f"{x['pct']['exact_at_rq2']:.1f}" for x in r["runs"]) + " |")


if __name__ == "__main__":
    main()
