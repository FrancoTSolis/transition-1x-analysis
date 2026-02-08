#!/usr/bin/env python
"""
Compare energy curves from DFT, ground truth (wB97x/6-31G(d)), and ANI-2x
side by side.  Ground truth is placed in the middle panel.

For each panel the energy barrier (E_TS - E_reactant) is annotated, where
the TS frame is defined as the frame with the highest ground-truth energy
among the selected frames.

On the DFT and ANI panels the ground-truth curve is overlaid and the
per-frame energy difference (method - GT) is highlighted with red hatched
shading.

Uses the same frame-selection logic as the animation scripts:
    [0, -last_frames, ..., -1, last_frames+1]

Usage:
    python scripts/compare_energy_methods.py \
        --xyz_dir trajectory_analysis_ani/C2H2N2O_rxn2091/xyz_frames \
        --ani_energies trajectory_analysis_ani/C2H2N2O_rxn2091/calculated_energies_ani.txt \
        --dft_cache trajectory_data.json \
        --last_frames 8

Requires: numpy, matplotlib (no torch / torchani / atkpython)
"""
import os
import re
import json
import argparse

import numpy as np
import matplotlib.pyplot as plt

HARTREE_TO_EV = 27.211386245988


# ---------------------------------------------------------------------------
# Frame selection (same logic as animation scripts)
# ---------------------------------------------------------------------------

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

        if 0 not in indices:
            extra_front.append(0)

        product_idx = last_k_frames + 1
        if product_idx < total_frames and product_idx not in indices:
            extra_back.append(product_idx)

        indices = extra_front + indices + extra_back

    return indices


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]


def read_gt_energies_from_xyz(xyz_dir):
    """Read ALL ground-truth energies (eV) from XYZ comment lines."""
    import glob
    frames = sorted(
        glob.glob(os.path.join(xyz_dir, "frame_*.xyz")),
        key=natural_sort_key,
    )
    if not frames:
        raise FileNotFoundError(f"No frame_*.xyz files found in {xyz_dir}")

    energies = []
    for fpath in frames:
        with open(fpath, 'r') as f:
            _ = f.readline()
            comment = f.readline()
        match = re.search(r'Energy=([-\d.eE+]+)', comment)
        if match:
            energies.append(float(match.group(1)))
        else:
            raise ValueError(
                f"No Energy= field in comment line of {fpath}."
            )
    return np.array(energies)


def read_energies_file(filepath, unit='hartree'):
    """Read energies from a plain-text file (one value per line)."""
    raw = np.loadtxt(filepath)
    if unit == 'hartree':
        return raw * HARTREE_TO_EV
    return raw


def read_dft_cache(filepath):
    """Read DFT energies from the JSON cache (animate_reaction.py format)."""
    with open(filepath, 'r') as f:
        data = json.load(f)

    entries = []
    for d in data:
        fname = os.path.basename(d['frame_file'])
        match = re.search(r'frame_(\d+)', fname)
        if match:
            idx = int(match.group(1))
            entries.append((idx, d['energy_ev']))

    entries.sort(key=lambda x: x[0])
    energies = np.array([e for _, e in entries])
    indices = [i for i, _ in entries]
    print(f"  {len(energies)} frames loaded (indices {indices[0]}\u2013{indices[-1]}).")
    return energies


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _annotate_barrier(ax, x_positions, energies_rel, ts_pos, color):
    """Draw a double-headed arrow for the energy barrier and label it."""
    barrier = energies_rel[ts_pos]
    x = x_positions[ts_pos]

    ax.annotate(
        '', xy=(x, barrier), xytext=(x, 0),
        arrowprops=dict(arrowstyle='<->', color=color, lw=1.8,
                        shrinkA=0, shrinkB=0),
    )
    ax.text(
        x + 0.3, barrier / 2,
        f' \u0394E\u2021 = {barrier:.3f} eV',
        fontsize=9, fontweight='bold', color=color,
        ha='left', va='center',
        bbox=dict(facecolor='white', edgecolor=color, alpha=0.85, pad=2),
    )
    ax.plot(x, barrier, '*', color=color, markersize=14,
            markeredgecolor='k', markeredgewidth=0.5, zorder=5)


def _setup_xaxis(ax, x_positions, frame_indices):
    """Set x-ticks to show actual frame indices."""
    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(i) for i in frame_indices], fontsize=7, rotation=45)
    ax.set_xlabel('Frame Index')


def _plot_gt_panel(ax, x_pos, gt_rel, ts_pos, frame_indices):
    """Plot the ground-truth panel (middle)."""
    ax.plot(x_pos, gt_rel, 'o-', color='#2ca02c', lw=1.8, markersize=5,
            label='GT: \u03c9B97x/6-31G(d)')
    ax.fill_between(x_pos, 0, gt_rel, alpha=0.10, color='#2ca02c')

    _annotate_barrier(ax, x_pos, gt_rel, ts_pos, '#2ca02c')

    ax.axhline(0, color='grey', lw=0.5, ls='--')
    ax.axvline(x_pos[ts_pos], color='grey', lw=0.5, ls=':')

    ax.set_title('Ground Truth \u2014 \u03c9B97x/6-31G(d)',
                 fontsize=11, fontweight='bold')
    ax.set_ylabel('Relative Energy (eV)')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    _setup_xaxis(ax, x_pos, frame_indices)


def _plot_method_panel(ax, x_pos, method_rel, gt_rel, ts_pos,
                       frame_indices, method_name, color):
    """Plot a DFT or ANI panel with error shading vs GT."""
    # GT reference (dashed)
    ax.plot(x_pos, gt_rel, 'o--', color='#2ca02c', lw=1.2, markersize=4,
            alpha=0.7, label='GT: \u03c9B97x/6-31G(d)')
    # Method curve
    ax.plot(x_pos, method_rel, 'o-', color=color, lw=1.8, markersize=5,
            label=method_name)

    # Red hatched fill for the difference
    ax.fill_between(
        x_pos, gt_rel, method_rel,
        facecolor='none', edgecolor='red', hatch='///',
        linewidth=0.0, alpha=0.6, label='Error vs GT',
    )
    ax.fill_between(x_pos, gt_rel, method_rel, color='red', alpha=0.08)

    _annotate_barrier(ax, x_pos, method_rel, ts_pos, color)

    ax.axhline(0, color='grey', lw=0.5, ls='--')
    ax.axvline(x_pos[ts_pos], color='grey', lw=0.5, ls=':')

    ax.set_title(method_name, fontsize=11, fontweight='bold')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    _setup_xaxis(ax, x_pos, frame_indices)


def _plot_empty_panel(ax, label):
    """Placeholder when data is unavailable."""
    ax.text(0.5, 0.5, f'{label}\n(no data)',
            ha='center', va='center', fontsize=14, color='grey',
            transform=ax.transAxes)
    ax.set_title(label, fontsize=11, fontweight='bold', color='grey')
    ax.set_xlabel('Frame Index')
    ax.grid(True, alpha=0.3)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compare DFT / Ground Truth / ANI-2x energy curves '
                    'side by side.',
    )
    parser.add_argument(
        '--xyz_dir', type=str, required=True,
        help='Directory containing frame_*.xyz files with GT energies',
    )
    parser.add_argument(
        '--ani_energies', type=str, default=None,
        help='Path to calculated_energies_ani.txt (Hartree)',
    )
    parser.add_argument(
        '--dft_energies', type=str, default=None,
        help='Path to calculated_energies.txt (Hartree, one per line)',
    )
    parser.add_argument(
        '--dft_cache', type=str, default=None,
        help='Path to trajectory_data.json (JSON cache, eV)',
    )
    parser.add_argument(
        '--last_frames', type=int, default=8,
        help='Number of last frames to include (default: 8)',
    )
    parser.add_argument(
        '--include_reactant_and_product', action=argparse.BooleanOptionalAction,
        default=True,
        help='Include frame 0 (reactant) and frame last_frames+1 (product). '
             'Default: on. Use --no-include_reactant_and_product to disable.',
    )
    parser.add_argument(
        '--output', type=str, default='energy_comparison_3way.png',
        help='Output image filename (default: energy_comparison_3way.png)',
    )
    args = parser.parse_args()

    # --- Load ALL energies ---
    print("Reading ground-truth energies from XYZ comment lines...")
    gt_all = read_gt_energies_from_xyz(args.xyz_dir)
    total_frames = len(gt_all)
    print(f"  {total_frames} total frames.")

    ani_all = None
    if args.ani_energies and os.path.exists(args.ani_energies):
        print(f"Reading ANI-2x energies from {args.ani_energies} (Hartree \u2192 eV)...")
        ani_all = read_energies_file(args.ani_energies, unit='hartree')
        print(f"  {len(ani_all)} frames loaded.")

    dft_all = None
    if args.dft_cache and os.path.exists(args.dft_cache):
        print(f"Reading DFT energies from {args.dft_cache} (JSON cache, eV)...")
        dft_all = read_dft_cache(args.dft_cache)
    elif args.dft_energies and os.path.exists(args.dft_energies):
        print(f"Reading DFT energies from {args.dft_energies} (Hartree \u2192 eV)...")
        dft_all = read_energies_file(args.dft_energies, unit='hartree')
        print(f"  {len(dft_all)} frames loaded.")

    # --- Select frames (same logic as animation scripts) ---
    indices = select_frame_indices(
        total_frames, args.last_frames, args.include_reactant_and_product,
    )
    print(f"\nSelected {len(indices)} frames: {indices}")

    # Subset energies
    gt_ev = gt_all[indices]

    ani_ev = None
    if ani_all is not None:
        ani_ev = ani_all[indices]

    dft_ev = None
    if dft_all is not None:
        dft_ev = dft_all[indices]

    # --- Relative energies (relative to first selected frame = reactant) ---
    gt_rel = gt_ev - gt_ev[0]

    ani_rel = None
    if ani_ev is not None:
        ani_rel = ani_ev - ani_ev[0]

    dft_rel = None
    if dft_ev is not None:
        dft_rel = dft_ev - dft_ev[0]

    # --- TS = highest-energy frame in selected GT ---
    ts_pos = int(np.argmax(gt_rel))     # position in the selected list
    ts_frame = indices[ts_pos]          # actual frame index
    barrier_gt = gt_rel[ts_pos]
    print(f"Ground-truth TS at frame {ts_frame} (position {ts_pos}), "
          f"barrier = {barrier_gt:.4f} eV")

    if ani_rel is not None:
        print(f"ANI-2x    barrier at same frame: {ani_rel[ts_pos]:.4f} eV "
              f"(error: {ani_rel[ts_pos] - barrier_gt:+.4f} eV)")

    if dft_rel is not None:
        print(f"DFT (PBE) barrier at same frame: {dft_rel[ts_pos]:.4f} eV "
              f"(error: {dft_rel[ts_pos] - barrier_gt:+.4f} eV)")

    # --- Plot ---
    x_pos = np.arange(len(indices))

    fig, (ax_dft, ax_gt, ax_ani) = plt.subplots(
        1, 3, figsize=(18, 6), sharey=True,
    )

    # Left: DFT
    if dft_rel is not None:
        _plot_method_panel(ax_dft, x_pos, dft_rel, gt_rel, ts_pos,
                           indices, 'DFT \u2014 GGA.PBE (QuantumATK)', '#1f77b4')
    else:
        _plot_empty_panel(ax_dft, 'DFT \u2014 GGA.PBE')

    # Middle: Ground Truth
    _plot_gt_panel(ax_gt, x_pos, gt_rel, ts_pos, indices)

    # Right: ANI
    if ani_rel is not None:
        _plot_method_panel(ax_ani, x_pos, ani_rel, gt_rel, ts_pos,
                           indices, 'ANI-2x', '#ff7f0e')
    else:
        _plot_empty_panel(ax_ani, 'ANI-2x')

    fig.suptitle(
        f'Energy Profile Comparison \u2014 {len(indices)} Selected Frames',
        fontsize=13, fontweight='bold', y=1.01,
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {args.output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
