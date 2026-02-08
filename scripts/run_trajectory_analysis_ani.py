#!/usr/bin/env python
"""
Run transition state trajectory analysis using the ANI-2x neural network potential.

This is the ANI equivalent of run_trajectory_analysis.py (which uses QuantumATK).
Unlike the ATK version, this script:
  - Runs with standard Python (not atkpython)
  - Computes all frames in a single batch (much faster)
  - Does not require subprocess calls

Requires: torch, torchani, h5py, numpy, matplotlib

Usage:
    python scripts/run_trajectory_analysis_ani.py [options]

Examples:
    # Single reaction (default)
    python scripts/run_trajectory_analysis_ani.py

    # Specific reaction
    python scripts/run_trajectory_analysis_ani.py --formula C2H2N2O --rxn_id rxn2091

    # Random sample of 5 reactions
    python scripts/run_trajectory_analysis_ani.py --formula "" --rxn_id "" --num_trajectories 5

    # All reactions
    python scripts/run_trajectory_analysis_ani.py --formula "" --rxn_id "" --all

    # Use GPU
    python scripts/run_trajectory_analysis_ani.py --device cuda
"""
import os
import sys
import json
import time
import random
import argparse

import numpy as np
import matplotlib.pyplot as plt
import torch
import torchani

# Add current directory to path to import extract_and_visualize
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import extract_and_visualize as ev


def load_ani_model(device="cpu"):
    """Load and configure the ANI-2x model."""
    print(f"Loading ANI-2x model on {device}...")
    model = torchani.models.ANI2x(periodic_table_index=True).to(device)
    model = model.float()
    for param in model.parameters():
        param.requires_grad = False
    print("ANI-2x model loaded successfully.")
    return model


def compute_energies_ani(model, atomic_numbers, positions, device="cpu"):
    """
    Compute energies for all frames using ANI-2x in a single batch.

    Args:
        model: ANI-2x model
        atomic_numbers: numpy array of shape [natoms] (atomic numbers)
        positions: numpy array of shape [nframes, natoms, 3] (Angstroms)
        device: torch device string

    Returns:
        numpy array of energies in Hartree, shape [nframes]
    """
    nframes = positions.shape[0]

    # Prepare tensors
    species = torch.tensor(atomic_numbers, dtype=torch.long, device=device)
    species_batch = species.unsqueeze(0).expand(nframes, -1)  # [nframes, natoms]

    pos_tensor = torch.tensor(positions, dtype=torch.float32, device=device)  # [nframes, natoms, 3]

    # Compute energies in batch
    with torch.no_grad():
        result = model((species_batch, pos_tensor))
        energies = result.energies.cpu().numpy()  # [nframes], in Hartree

    # Check for NaN/Inf
    bad_mask = np.isnan(energies) | np.isinf(energies)
    if bad_mask.any():
        bad_indices = np.where(bad_mask)[0]
        print(f"  Warning: NaN/Inf energy at frames: {bad_indices.tolist()}")
        energies[bad_mask] = np.nan

    return energies


def process_trajectory(f, formula, rxn_id, args, model, device):
    """
    Process a single trajectory defined by formula and rxn_id.
    """
    print(f"\nProcessing Reaction: {formula} - {rxn_id}")

    # Create output directory for this specific reaction
    reaction_dir = os.path.join(args.output_dir, f"{formula}_{rxn_id}")
    os.makedirs(reaction_dir, exist_ok=True)
    xyz_dir = os.path.join(reaction_dir, "xyz_frames")
    os.makedirs(xyz_dir, exist_ok=True)

    try:
        grp, _, _ = ev.get_reaction_data(f, formula, rxn_id)
        z, pos, original_energies = ev.extract_trajectory(grp)
    except Exception as e:
        print(f"Error loading data for {formula}-{rxn_id}: {e}")
        return

    nframes = pos.shape[0]

    # Save XYZ files for each frame
    for i in range(nframes):
        xyz_file = os.path.join(xyz_dir, f"frame_{i:03d}.xyz")
        ev.save_xyz(xyz_file, z, pos[i], original_energies[i])

    if args.dry_run:
        print(f"  Dry run: Skipping ANI-2x calculations for {nframes} frames.")
        calculated_energies = original_energies.copy()
    else:
        # Compute all energies in one batch — much faster than per-frame subprocess calls
        print(f"  Computing ANI-2x energies for {nframes} frames in batch...", end="", flush=True)
        t0 = time.time()
        calculated_energies = compute_energies_ani(model, z, pos, device)
        dt = time.time() - t0
        print(f" Done ({dt:.2f}s)")

        # Print per-frame summary
        for i in range(nframes):
            if not np.isnan(calculated_energies[i]):
                print(f"    Frame {i:3d}: E = {calculated_energies[i]:.6f} Ha")
            else:
                print(f"    Frame {i:3d}: FAILED")

    # --- Plot Results ---
    plt.figure(figsize=(10, 6))

    orig_rel = original_energies - original_energies[0]

    valid_calc = np.array(calculated_energies)
    if not np.all(np.isnan(valid_calc)):
        first_valid_idx = np.where(~np.isnan(valid_calc))[0][0]
        calc_rel = valid_calc - valid_calc[first_valid_idx]
        # Align to original at the first valid point
        calc_rel = calc_rel - calc_rel[first_valid_idx] + orig_rel[first_valid_idx]
    else:
        calc_rel = valid_calc

    plt.plot(range(len(original_energies)), orig_rel, 'o-', label='Original (wB97x/6-31G(d))', alpha=0.7)
    plt.plot(range(len(calculated_energies)), calc_rel, 's--', label='Calculated (ANI-2x)', alpha=0.7)

    plt.xlabel('Frame Index')
    plt.ylabel('Relative Energy (Hartree)')
    plt.title(f'Energy Profile Comparison: {formula} - {rxn_id}')
    plt.legend()
    plt.grid(True)

    plot_path = os.path.join(reaction_dir, "energy_comparison_ani.png")
    plt.savefig(plot_path)
    plt.close()

    # Save raw data
    np.savetxt(os.path.join(reaction_dir, "calculated_energies_ani.txt"), calculated_energies)

    # Save a comparison summary as JSON
    summary = {
        "formula": formula,
        "rxn_id": rxn_id,
        "nframes": int(nframes),
        "method": "ANI-2x",
        "original_energies_hartree": original_energies.tolist(),
        "ani_energies_hartree": calculated_energies.tolist(),
    }
    with open(os.path.join(reaction_dir, "summary_ani.json"), 'w') as fp:
        json.dump(summary, fp, indent=2)

    print(f"  Analysis saved to {reaction_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Run trajectory analysis with ANI-2x neural network potential"
    )
    parser.add_argument(
        "--h5_path", type=str,
        default="Transition-State-Generation-Flow/dataset/transition1x/data/Transition1x.h5",
        help="Path to h5 file",
    )
    parser.add_argument(
        "--output_dir", type=str, default="trajectory_analysis_ani",
        help="Output directory",
    )
    parser.add_argument(
        "--formula", type=str, default="C2H2N2O",
        help="Chemical formula to select",
    )
    parser.add_argument(
        "--rxn_id", type=str, default="rxn2091",
        help="Reaction ID to select",
    )
    parser.add_argument(
        "--num_trajectories", type=int, default=1,
        help="Number of random trajectories to sample if formula/rxn_id not specified",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Run all trajectories in the dataset (overrides num_trajectories)",
    )
    parser.add_argument(
        "--dry_run", action="store_true",
        help="Don't actually run ANI, just extract and setup",
    )
    parser.add_argument(
        "--device", type=str, default="cpu", choices=["cpu", "cuda"],
        help="Device to run ANI on (default: cpu)",
    )

    args = parser.parse_args()

    # Load dataset
    try:
        f = ev.load_data(args.h5_path)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # Load ANI model once (skip in dry_run)
    device = torch.device(args.device)
    if not args.dry_run:
        model = load_ani_model(device)
    else:
        model = None

    # Determine which trajectories to process
    trajectories_to_process = []

    if args.formula and args.rxn_id:
        # Specific trajectory requested
        trajectories_to_process.append((args.formula, args.rxn_id))
    else:
        # Get all available trajectories
        all_reactions = ev.get_all_reactions(f)
        total_reactions = len(all_reactions)
        print(f"Found {total_reactions} total reactions in dataset.")

        if args.all:
            trajectories_to_process = all_reactions
        else:
            num_to_sample = min(args.num_trajectories, total_reactions)
            trajectories_to_process = random.sample(all_reactions, num_to_sample)
            print(f"Randomly selected {num_to_sample} trajectories.")

    # Process each trajectory
    print(f"Processing {len(trajectories_to_process)} trajectories...")
    for i, (formula, rxn_id) in enumerate(trajectories_to_process):
        print(f"[{i+1}/{len(trajectories_to_process)}]", end=" ")
        process_trajectory(f, formula, rxn_id, args, model, device)


if __name__ == "__main__":
    main()
