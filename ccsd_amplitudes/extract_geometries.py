#!/usr/bin/env python3
"""Extract reactant / TS / product geometries from Transition1x HDF5
and generate Q-Chem CCSD input files with PRINT_QIS for full amplitude dump.

Usage:
    python extract_geometries.py                         # all reactions
    python extract_geometries.py --max-atoms 14          # filter by size
    python extract_geometries.py --formula C2H2N2O       # single formula
    python extract_geometries.py --dry-run               # preview without writing
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

PERIODIC_TABLE = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 35: "Br", 53: "I",
}

ELECTRON_COUNT = {1: 1, 6: 6, 7: 7, 8: 8, 9: 9, 16: 16, 17: 17, 15: 15, 35: 35, 53: 53}
BF_STO3G = {1: 1, 6: 5, 7: 5, 8: 5, 9: 5, 16: 9, 17: 9, 15: 9, 35: 14, 53: 19}

QCHEM_TEMPLATE = """\
$molecule
0 1
{coords}
$end

$rem
    METHOD          ccsd
    BASIS           STO-3G
    UNRESTRICTED    TRUE
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   DIIS_GDM
    MAX_SCF_CYCLES  200
    THRESH          15
    SCF_CONVERGENCE 10
    CC_CONVERGENCE  7
    N_FROZEN_CORE   FC
    GUI             2
    PRINT_QIS       TRUE
$end
"""


def format_coords(atomic_numbers: np.ndarray, positions: np.ndarray) -> str:
    lines = []
    for z, pos in zip(atomic_numbers, positions):
        sym = PERIODIC_TABLE.get(int(z), "X")
        lines.append(f"{sym:2s}  {pos[0]:14.8f}  {pos[1]:14.8f}  {pos[2]:14.8f}")
    return "\n".join(lines)


def write_xyz(path: Path, atomic_numbers: np.ndarray, positions: np.ndarray,
              comment: str = "") -> None:
    with open(path, "w") as fh:
        fh.write(f"{len(atomic_numbers)}\n")
        fh.write(f"{comment}\n")
        for z, pos in zip(atomic_numbers, positions):
            sym = PERIODIC_TABLE.get(int(z), "X")
            fh.write(f"{sym} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n")


def write_qchem_input(path: Path, atomic_numbers: np.ndarray,
                      positions: np.ndarray) -> None:
    coords = format_coords(atomic_numbers, positions)
    content = QCHEM_TEMPLATE.format(coords=coords)
    path.write_text(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5-path", default="../data/transition1x.h5",
                        help="Path to Transition1x HDF5 file")
    parser.add_argument("--output-dir", default="jobs",
                        help="Output directory for Q-Chem inputs and XYZ files")
    parser.add_argument("--formula", default=None,
                        help="Process only this chemical formula")
    parser.add_argument("--rxn-id", default=None,
                        help="Process only this reaction ID (requires --formula)")
    parser.add_argument("--max-atoms", type=int, default=None,
                        help="Skip molecules with more atoms than this")
    parser.add_argument("--min-atoms", type=int, default=None,
                        help="Skip molecules with fewer atoms than this")
    parser.add_argument("--max-reactions", type=int, default=None,
                        help="Limit total number of reactions to process")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print summary without writing files")
    parser.add_argument("--split", default="data",
                        choices=["data", "train", "test", "val"],
                        help="Which HDF5 split to use (default: data = union)")
    args = parser.parse_args()

    h5_path = Path(args.h5_path)
    if not h5_path.exists():
        print(f"ERROR: HDF5 file not found: {h5_path}")
        sys.exit(1)

    f = h5py.File(str(h5_path), "r")
    data = f[args.split]

    formulas = [args.formula] if args.formula else sorted(data.keys())
    output_dir = Path(args.output_dir)

    total = 0
    skipped_size = 0
    manifest = []

    for formula in formulas:
        if formula not in data:
            print(f"WARNING: formula {formula} not found in split '{args.split}'")
            continue

        rxn_ids = [args.rxn_id] if args.rxn_id else sorted(data[formula].keys())

        for rxn_id in rxn_ids:
            if rxn_id not in data[formula]:
                print(f"WARNING: {formula}/{rxn_id} not found")
                continue

            grp = data[formula][rxn_id]
            z = grp["atomic_numbers"][:]
            natoms = len(z)

            if args.max_atoms and natoms > args.max_atoms:
                skipped_size += 1
                continue
            if args.min_atoms and natoms < args.min_atoms:
                skipped_size += 1
                continue

            n_elec = sum(ELECTRON_COUNT.get(int(a), 0) for a in z)
            n_bf = sum(BF_STO3G.get(int(a), 0) for a in z)
            n_occ = n_elec // 2
            n_virt = n_bf - n_occ

            entry = {
                "formula": formula,
                "rxn_id": rxn_id,
                "natoms": natoms,
                "n_electrons": n_elec,
                "n_basis_sto3g": n_bf,
                "n_occ": n_occ,
                "n_virt": n_virt,
            }

            for label in ["reactant", "transition_state", "product"]:
                sub = grp[label]
                pos = sub["positions"][0]
                energy = float(sub["wB97x_6-31G(d).energy"][0])
                entry[f"{label}_energy_eV"] = energy

                if not args.dry_run:
                    tag = {"reactant": "R", "transition_state": "TS",
                           "product": "P"}[label]
                    job_name = f"{formula}_{rxn_id}_{tag}"
                    job_dir = output_dir / job_name
                    job_dir.mkdir(parents=True, exist_ok=True)

                    write_xyz(job_dir / f"{job_name}.xyz", z, pos,
                              comment=f"{formula} {rxn_id} {label} E={energy:.6f} eV")
                    write_qchem_input(job_dir / f"{job_name}.in", z, pos)

            manifest.append(entry)
            total += 1

            if args.max_reactions and total >= args.max_reactions:
                break
        if args.max_reactions and total >= args.max_reactions:
            break

    f.close()

    print(f"\nProcessed {total} reactions ({skipped_size} skipped by size filter)")
    print(f"Total Q-Chem jobs to run: {total * 3} (3 per reaction: R, TS, P)")

    if manifest:
        natoms_list = [e["natoms"] for e in manifest]
        nbf_list = [e["n_basis_sto3g"] for e in manifest]
        print(f"Atom count range: {min(natoms_list)}–{max(natoms_list)}")
        print(f"Basis function range (STO-3G): {min(nbf_list)}–{max(nbf_list)}")

        # Estimate disk usage for amplitude files
        total_t1_elements = 0
        total_t2_elements = 0
        for e in manifest:
            nocc, nvirt = e["n_occ"], e["n_virt"]
            total_t1_elements += nocc * nvirt
            total_t2_elements += nocc * nocc * nvirt * nvirt
        # 3 species per reaction, ~25 bytes per element in text format
        bytes_per_element = 25
        per_rxn_t1 = sum(e["n_occ"] * e["n_virt"] for e in manifest) / len(manifest) * bytes_per_element * 3
        per_rxn_t2 = sum(e["n_occ"]**2 * e["n_virt"]**2 for e in manifest) / len(manifest) * bytes_per_element * 3
        total_amp_bytes = (total_t1_elements + total_t2_elements) * bytes_per_element * 3
        # Scratch files: ~15x the amplitude file size based on the 14-atom reference
        total_scratch_bytes = total_amp_bytes * 15

        print(f"\nDisk space estimates (amplitude .dat files only):")
        print(f"  T1 total: {total_t1_elements * bytes_per_element * 3 / 1e9:.2f} GB")
        print(f"  T2 total: {total_t2_elements * bytes_per_element * 3 / 1e9:.2f} GB")
        print(f"  Combined amplitude files: {total_amp_bytes / 1e9:.2f} GB")
        print(f"  Estimated total scratch (during run): {total_scratch_bytes / 1e9:.1f} GB")
        print(f"  Average per reaction (3 species): {(per_rxn_t1 + per_rxn_t2) / 1e6:.1f} MB amp files")

    if not args.dry_run:
        manifest_path = output_dir / "manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        print(f"\nManifest written to {manifest_path}")
        print(f"Job directories created under {output_dir}/")
    else:
        print("\n(Dry run — no files written)")


if __name__ == "__main__":
    main()
