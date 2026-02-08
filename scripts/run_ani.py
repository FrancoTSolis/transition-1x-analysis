#!/usr/bin/env python
"""
Calculate energy of a molecule from an XYZ file using the ANI-2x neural network potential.

This is the ANI equivalent of run_dft.py (which uses QuantumATK).
Unlike run_dft.py, this script runs with standard Python (not atkpython).

Requires: torch, torchani

Usage:
    python scripts/run_ani.py <xyz_file> [--device cpu|cuda]

Supported elements (ANI-2x): H, C, N, O, F, S, Cl
"""
import sys
import os
import json
import argparse

import torch
import torchani

# Element symbol to atomic number mapping (ANI-2x supported elements only)
SYMBOL_TO_Z = {
    'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16, 'Cl': 17,
}
Z_TO_SYMBOL = {v: k for k, v in SYMBOL_TO_Z.items()}

ANI2X_SUPPORTED = set(SYMBOL_TO_Z.keys())

# Conversion factor
HARTREE_TO_EV = 27.211386245988


def read_xyz(filename):
    """
    Read an XYZ file and return atomic numbers and positions.

    Returns:
        atomic_numbers: list of int
        positions: list of [x, y, z] floats (Angstroms)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    try:
        num_atoms = int(lines[0].strip())
    except (ValueError, IndexError):
        raise ValueError("Invalid XYZ file format: first line must be number of atoms")

    atomic_numbers = []
    positions = []

    for i, line in enumerate(lines[2:]):
        if i >= num_atoms:
            break
        parts = line.split()
        if len(parts) < 4:
            continue

        symbol = parts[0]
        if symbol not in SYMBOL_TO_Z:
            raise ValueError(
                f"Element '{symbol}' is not supported by ANI-2x. "
                f"Supported elements: {sorted(ANI2X_SUPPORTED)}"
            )

        atomic_numbers.append(SYMBOL_TO_Z[symbol])
        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if len(atomic_numbers) != num_atoms:
        print(f"Warning: Expected {num_atoms} atoms, found {len(atomic_numbers)}")

    return atomic_numbers, positions


def compute_energy(model, atomic_numbers, positions, device):
    """
    Compute energy for a single molecule using ANI-2x.

    Args:
        model: loaded ANI-2x model
        atomic_numbers: list of int
        positions: list of [x, y, z]
        device: torch.device

    Returns:
        energy in Hartree (float)
    """
    species = torch.tensor([atomic_numbers], dtype=torch.long, device=device)   # [1, natoms]
    coords = torch.tensor([positions], dtype=torch.float32, device=device)      # [1, natoms, 3]

    with torch.no_grad():
        result = model((species, coords))
        energy_hartree = result.energies.item()

    return energy_hartree


def main():
    parser = argparse.ArgumentParser(
        description="Calculate energy using ANI-2x neural network potential"
    )
    parser.add_argument("xyz_file", type=str, help="Path to XYZ file")
    parser.add_argument(
        "--device", type=str, default="cpu", choices=["cpu", "cuda"],
        help="Device to run on (default: cpu)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.xyz_file):
        print(f"Error: File {args.xyz_file} not found.")
        sys.exit(1)

    print(f"Running ANI-2x for {args.xyz_file}...")

    # 1. Read XYZ
    try:
        atomic_numbers, positions = read_xyz(args.xyz_file)
    except Exception as e:
        print(f"Error reading XYZ file: {e}")
        sys.exit(1)

    # 2. Load ANI-2x model
    device = torch.device(args.device)
    print(f"Loading ANI-2x model on {device}...")
    model = torchani.models.ANI2x(periodic_table_index=True).to(device)
    model = model.float()

    # 3. Compute energy
    energy_hartree = compute_energy(model, atomic_numbers, positions, device)
    energy_ev = energy_hartree * HARTREE_TO_EV

    print(f"Energy: {energy_hartree:.8f} Hartree  ({energy_ev:.6f} eV)")

    # 4. Output results (compatible format with run_dft.py)
    output_data = {
        "energy_hartree": energy_hartree,
        "energy_ev": energy_ev,
        "method": "ANI-2x",
        "atoms": [
            {"symbol": Z_TO_SYMBOL[z], "atomic_number": z}
            for z in atomic_numbers
        ],
    }

    print("JSON_OUTPUT_START")
    print(json.dumps(output_data))
    print("JSON_OUTPUT_END")

    # Compatible single-line output (mirrors DFT_ENERGY_HARTREE from run_dft.py)
    print(f"ANI_ENERGY_HARTREE: {energy_hartree}")


if __name__ == "__main__":
    main()
