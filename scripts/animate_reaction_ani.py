#!/usr/bin/env python
"""
Animate a reaction trajectory using ANI-2x energies.

This is the ANI equivalent of animate_reaction.py (which uses QuantumATK DFT).
Unlike the DFT version:
  - Uses ANI-2x for energy calculation (no atkpython needed)
  - Computes all frame energies in a single batch (much faster)
  - Colors atoms by element (ANI does not provide partial charges)

Usage:
    python scripts/animate_reaction_ani.py --frames_dir trajectory_analysis_ani/C2H2N2O_rxn2091/xyz_frames --last_frames 8

Requires: torch, torchani, numpy, matplotlib, scipy
"""
import os
import re
import json
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.interpolate import make_interp_spline

import torch
import torchani

# ANI-2x supported elements
SYMBOL_TO_Z = {
    'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16, 'Cl': 17,
}

# Element colors (CPK convention)
ELEMENT_COLORS = {
    'H': '#FFFFFF',   # white
    'C': '#909090',   # grey
    'N': '#3050F8',   # blue
    'O': '#FF0D0D',   # red
    'F': '#90E050',   # green
    'S': '#FFFF30',   # yellow
    'Cl': '#1FF01F',  # green
}

# Element sizes (relative)
ELEMENT_SIZES = {
    'H': 50, 'C': 100, 'N': 100, 'O': 100, 'F': 90, 'S': 120, 'Cl': 110,
}

# Conversion
HARTREE_TO_EV = 27.211386245988


def natural_sort_key(s):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]


def read_xyz_coords(filename):
    """Read an XYZ file and return symbols and coordinates."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip())
    coords = []
    symbols = []
    for line in lines[2:2 + num_atoms]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


def compute_ani_energies(frames_data, device="cpu"):
    """
    Compute ANI-2x energies for all frames in a single batch.

    Args:
        frames_data: list of dicts with 'symbols' and 'coordinates'
        device: torch device string

    Returns:
        list of energies in eV
    """
    print(f"Loading ANI-2x model on {device}...")
    model = torchani.models.ANI2x(periodic_table_index=True).to(device)
    model = model.float()

    nframes = len(frames_data)

    # Convert symbols to atomic numbers
    atomic_numbers = [SYMBOL_TO_Z[s] for s in frames_data[0]['symbols']]
    species = torch.tensor(atomic_numbers, dtype=torch.long, device=device)
    species_batch = species.unsqueeze(0).expand(nframes, -1)  # [nframes, natoms]

    # Stack all coordinates
    all_coords = np.array([d['coordinates'] for d in frames_data])  # [nframes, natoms, 3]
    pos_tensor = torch.tensor(all_coords, dtype=torch.float32, device=device)

    # Batch compute
    print(f"Computing ANI-2x energies for {nframes} frames...", end="", flush=True)
    with torch.no_grad():
        result = model((species_batch, pos_tensor))
        energies_hartree = result.energies.cpu().numpy()

    energies_ev = energies_hartree * HARTREE_TO_EV
    print(" Done.")

    return energies_ev.tolist(), energies_hartree.tolist()


def select_frame_indices(total_frames, last_k_frames, include_reactant_and_product):
    """
    Build the ordered list of frame indices to include.

    Default (include_reactant_and_product=True):
        [0, -last_k, -(last_k-1), ..., -1, last_k+1]
        i.e. reactant, last K frames, product

    Without the flag:
        just the last K frames (original behaviour).
    """
    if last_k_frames > 0:
        start_idx = max(0, total_frames - last_k_frames)
        indices = list(range(start_idx, total_frames))
    else:
        indices = list(range(total_frames))

    if include_reactant_and_product:
        extra_front = []
        extra_back = []

        # Frame 0 = reactant
        if 0 not in indices:
            extra_front.append(0)

        # Frame (last_k_frames + 1) = product
        product_idx = last_k_frames + 1
        if product_idx < total_frames and product_idx not in indices:
            extra_back.append(product_idx)

        indices = extra_front + indices + extra_back

    return indices


def load_frames(frame_dir, last_k_frames=10, include_reactant_and_product=True):
    """
    Load XYZ frames from a directory.

    Args:
        frame_dir: path to directory containing frame_*.xyz files
        last_k_frames: number of last frames to use
        include_reactant_and_product: if True, also include frame 0 (reactant)
            at the front and frame last_k_frames+1 (product) at the end

    Returns:
        list of dicts with 'frame_file', 'symbols', 'coordinates'
    """
    import glob
    all_frames = sorted(glob.glob(os.path.join(frame_dir, "frame_*.xyz")), key=natural_sort_key)
    if not all_frames:
        print(f"No frames found in {frame_dir}")
        return []

    total_frames = len(all_frames)
    indices = select_frame_indices(total_frames, last_k_frames, include_reactant_and_product)

    desc_parts = []
    if include_reactant_and_product:
        desc_parts.append(f"frame 0 (reactant)")
    desc_parts.append(f"last {min(last_k_frames, total_frames)} frames")
    if include_reactant_and_product:
        product_idx = last_k_frames + 1
        if product_idx < total_frames:
            desc_parts.append(f"frame {product_idx} (product)")
    print(f"Selected {len(indices)} frames: {', '.join(desc_parts)}")
    print(f"  Indices: {indices}")

    data = []
    for idx in indices:
        frame_path = all_frames[idx]
        symbols, coords = read_xyz_coords(frame_path)
        data.append({
            'frame_file': frame_path,
            'frame_index': idx,
            'symbols': symbols,
            'coordinates': coords.tolist(),
        })

    return data


def create_animation(data, output_file='reaction_animation_ani.gif'):
    """Create animation GIF from frame data with energies."""
    if not data:
        print("No data to animate.")
        return

    # Extract energy array
    energies = np.array([d['energy_ev'] for d in data])
    energies_rel = energies - energies[0]  # Relative to first frame

    num_points = len(energies_rel)
    reaction_progress = np.linspace(0, 1, num_points)

    # Smooth curve via spline interpolation
    reaction_progress_smooth = np.linspace(0, 1, 200)
    k_order = min(3, num_points - 1)  # cubic if enough points
    if k_order >= 1:
        spl = make_interp_spline(reaction_progress, energies_rel, k=k_order)
        energies_smooth = spl(reaction_progress_smooth)
    else:
        reaction_progress_smooth = reaction_progress
        energies_smooth = energies_rel

    # --- Setup figure ---
    fig = plt.figure(figsize=(12, 6))
    ax3d = fig.add_subplot(1, 2, 1, projection='3d')
    ax_energy = fig.add_subplot(1, 2, 2)

    # Energy profile (background)
    ax_energy.plot(reaction_progress_smooth, energies_smooth, 'b-', label='Energy Profile', alpha=0.7)
    ax_energy.plot(reaction_progress, energies_rel, 'b.', alpha=0.3)
    ax_energy.set_xlabel(r'Reaction Progress $\xi$')
    ax_energy.set_ylabel('Relative Energy (eV)')
    ax_energy.set_title('Reaction Energy Profile (ANI-2x)')
    energy_dot, = ax_energy.plot([], [], 'ro', markersize=10)

    # Pre-calculate 3D bounds
    all_coords = np.array([d['coordinates'] for d in data])
    min_coord = all_coords.min(axis=(0, 1))
    max_coord = all_coords.max(axis=(0, 1))
    center = (min_coord + max_coord) / 2
    span = (max_coord - min_coord).max() / 2 + 0.5  # small padding

    def init():
        energy_dot.set_data([], [])
        return energy_dot,

    def update(frame_idx):
        frame_data = data[frame_idx]
        coords = np.array(frame_data['coordinates'])
        symbols = frame_data['symbols']

        ax3d.clear()

        # Atom sizes and colors by element
        sizes = [ELEMENT_SIZES.get(s, 100) for s in symbols]
        colors = [ELEMENT_COLORS.get(s, '#808080') for s in symbols]

        ax3d.scatter(
            coords[:, 0], coords[:, 1], coords[:, 2],
            s=sizes, c=colors, depthshade=True, edgecolors='k', linewidths=0.5,
        )

        # Draw bonds (simple distance criterion)
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < 1.8:
                    ax3d.plot(
                        [coords[i, 0], coords[j, 0]],
                        [coords[i, 1], coords[j, 1]],
                        [coords[i, 2], coords[j, 2]],
                        'k-', linewidth=2, alpha=0.5,
                    )

        ax3d.set_xlim(center[0] - span, center[0] + span)
        ax3d.set_ylim(center[1] - span, center[1] + span)
        ax3d.set_zlim(center[2] - span, center[2] + span)
        ax3d.set_title(f'Frame {frame_idx} | E: {energies_rel[frame_idx]:.3f} eV')

        # Update energy dot
        current_progress = frame_idx / (len(data) - 1) if len(data) > 1 else 0
        energy_dot.set_data([current_progress], [energies_rel[frame_idx]])

        return ax3d, energy_dot

    # Build legend for element colors
    from matplotlib.lines import Line2D
    unique_symbols = sorted(set(data[0]['symbols']))
    legend_handles = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=ELEMENT_COLORS.get(s, '#808080'),
               markeredgecolor='k', markersize=8, label=s)
        for s in unique_symbols
    ]
    ax3d.legend(handles=legend_handles, loc='upper left', fontsize=8)

    ani = FuncAnimation(fig, update, frames=len(data), init_func=init, blit=False, interval=500)

    print(f"Saving animation to {output_file}...")
    try:
        ani.save(output_file, writer='pillow', fps=2)
        print(f"Done! Animation saved to {output_file}")
    except Exception as e:
        print(f"Error saving animation: {e}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Animate reaction trajectory with ANI-2x energies.'
    )
    parser.add_argument(
        '--frames_dir', type=str,
        default='trajectory_analysis_ani/C2H2N2O_rxn2091/xyz_frames',
        help='Directory containing XYZ frames',
    )
    parser.add_argument(
        '--last_frames', type=int, default=10,
        help='Number of last frames to process (default: 10)',
    )
    parser.add_argument(
        '--output', type=str, default='reaction_animation_ani.gif',
        help='Output GIF filename',
    )
    parser.add_argument(
        '--device', type=str, default='cpu', choices=['cpu', 'cuda'],
        help='Device to run ANI on (default: cpu)',
    )
    parser.add_argument(
        '--include_reactant_and_product', action=argparse.BooleanOptionalAction,
        default=True,
        help='Include frame 0 (reactant) and frame last_frames+1 (product) '
             'in the animation (default: on). Use --no-include_reactant_and_product to disable.',
    )

    args = parser.parse_args()

    # 1. Load frames
    data = load_frames(
        args.frames_dir,
        last_k_frames=args.last_frames,
        include_reactant_and_product=args.include_reactant_and_product,
    )
    if not data:
        return

    # 2. Compute ANI-2x energies for all frames
    energies_ev, energies_hartree = compute_ani_energies(data, device=args.device)

    # Attach energies to frame data
    for i, d in enumerate(data):
        d['energy_ev'] = energies_ev[i]
        d['energy_hartree'] = energies_hartree[i]

    # 3. Create animation
    create_animation(data, output_file=args.output)


if __name__ == "__main__":
    main()
