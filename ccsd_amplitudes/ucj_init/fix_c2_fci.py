"""Post-correct the C2 FCI reference in the small-system result files.

The first runs took E_FCI from pyscf's direct_spin1 with nroots=1, which for C2 (1.25 A, STO-3G, all electrons) lands
on an excited state: -74.642955 Eh. The ground state is -74.690586 Eh (nroots=6 and a Lanczos solve with
exact_model.DenseFCI agree). H4 and LiH were right. common.ref_energies now uses nroots=4.
This script rewrites e_fci (and the quantities derived from it) for C2 records; it is idempotent.

    python -m ucj_init.fix_c2_fci
"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / "pretrain" / "opt_true" / "results" / "ucj_init"
E_FCI_C2 = -74.690586


def pct(E, e_hf):
    return 100.0 * (e_hf - E) / (e_hf - E_FCI_C2)


def main():
    global E_HF_C2
    f = RES / "bch2_reproduce.json"
    rows = json.load(open(f))
    E_HF_C2 = next(r["e_hf"] for r in rows if r["name"].startswith("C2"))
    for r in rows:
        if r["name"].startswith("C2"):
            r["e_fci_wrong_root"] = r.get("e_fci_wrong_root", r["e_fci"])
            r["e_fci"] = E_FCI_C2
            r["fci_corr"] = E_FCI_C2 - r["e_hf"]
    json.dump(rows, open(f, "w"), indent=1)
    f = RES / "avijit_exact.json"
    if f.exists():
        rows = json.load(open(f))
        for r in rows:
            if r["name"].startswith("C2"):
                r["e_fci_wrong_root"] = r.get("e_fci_wrong_root", r["e_fci"])
                r["e_fci"] = E_FCI_C2
        json.dump(rows, open(f, "w"), indent=1)
    for name in ("ansatz_capacity.json", "ansatz_capacity_c2.json"):
        f = RES / name
        if not f.exists():
            continue
        rows = json.load(open(f))
        for r in rows:
            if not r["name"].startswith("C2"):
                continue
            e_hf = E_HF_C2
            r["e_hf"] = e_hf
            r["fci_corr_wrong_root"] = r.get("fci_corr_wrong_root", r["fci_corr"])
            r["fci_corr"] = E_FCI_C2 - e_hf
            for x in r["runs"]:
                x["pct"] = pct(x["E"], e_hf)
                x["pct0"] = pct(x["E0"], e_hf)
            r["best_pct"] = max(x["pct"] for x in r["runs"])
        json.dump(rows, open(f, "w"), indent=1)
    print("C2 FCI reference corrected to", E_FCI_C2)


if __name__ == "__main__":
    main()
