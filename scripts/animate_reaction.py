import os
import glob
import subprocess
import json
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]

def read_xyz_coords(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip())
    coords = []
    symbols = []
    for line in lines[2:2+num_atoms]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)

import argparse

def run_dft_on_frames(frame_dir, last_k_frames=10, cache_file='trajectory_data.json'):
    # Check for cache
    # Note: Cache might be stale if we change sampling parameters, 
    # but for simplicity we rely on user to clean it or we could make cache key depend on params.
    # Actually, if we just want last K frames, using a cached list of ALL frames (if we had it) would be fine,
    # but here the cache stores the RESULTS of the run. 
    # If the cache contains results for the frames we want, we could use them.
    # But simplifying: if cache exists, warn user or just overwrite?
    # The previous implementation loaded cache blindly. 
    # Let's force re-calculation if cache doesn't match roughly, 
    # OR better: just load cache, check if we have data for the specific frames we want.
    # For now, to ensure we get the right K frames, let's ignore cache for frame selection 
    # but maybe use it for lookup?
    # Simpler: If user asks for last K frames, we just pick them.
    
    frames = sorted(glob.glob(os.path.join(frame_dir, "frame_*.xyz")), key=natural_sort_key)
    if not frames:
        print(f"No frames found in {frame_dir}")
        return []
        
    # Select last K frames
    total_frames = len(frames)
    if last_k_frames > 0:
        start_idx = max(0, total_frames - last_k_frames)
        frames = frames[start_idx:]
        print(f"Selected last {len(frames)} frames (from index {start_idx} to {total_frames-1})")
    else:
        # If 0 or negative, maybe process all? Or just default behavior.
        # Let's assume argument is always provided with default 10.
        pass

    results = []
    # Load cache for lookup
    cache_data = {}
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r') as f:
                cached_list = json.load(f)
                # Convert list to dict keyed by filename
                for item in cached_list:
                    if 'frame_file' in item:
                        cache_data[item['frame_file']] = item
        except:
            print("Warning: Could not read cache file.")

    print(f"Processing {len(frames)} frames...")

    for i, frame in enumerate(frames):
        # Check cache first
        if frame in cache_data:
            print(f"Skipping {frame} (found in cache)")
            results.append(cache_data[frame])
            continue

        print(f"Processing frame {i+1}/{len(frames)}: {frame}")
        
        # Run DFT
        cmd = ["atkpython", "scripts/run_dft.py", frame]
        try:
            # Capture output
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse JSON from stdout
            # The script output contains other text, so we look for JSON block
            stdout = result.stdout
            json_start = stdout.find("JSON_OUTPUT_START")
            json_end = stdout.find("JSON_OUTPUT_END")
            
            if json_start != -1 and json_end != -1:
                json_str = stdout[json_start + len("JSON_OUTPUT_START"):json_end].strip()
                data = json.loads(json_str)
                
                # Add frame info
                symbols, coords = read_xyz_coords(frame)
                data['frame_file'] = frame
                data['frame_index'] = i # Relative index in this run
                data['coordinates'] = coords.tolist() # Store as list for JSON
                data['symbols'] = symbols # Redundant but safe
                
                results.append(data)
                
                # Update cache incrementally (optional, but good for safety)
                # For now just keep in memory and write at end or use existing logic
            else:
                print(f"Error parsing output for {frame}: JSON block not found.")
                print("Stdout:", stdout)

        except subprocess.CalledProcessError as e:
            print(f"Error running DFT on {frame}: {e}")
            print("Stderr:", e.stderr)
    
    # Save cache (merge with existing or overwrite? simpler to overwrite with current set 
    # BUT we might lose expensive calculations for other frames.
    # Better: Update the cache file with new results.
    for r in results:
        cache_data[r['frame_file']] = r
    
    with open(cache_file, 'w') as f:
        json.dump(list(cache_data.values()), f, indent=2)
    
    # Return only the results for the selected frames, in order
    ordered_results = []
    for frame in frames:
        if frame in cache_data:
            ordered_results.append(cache_data[frame])
            
    return ordered_results

def create_animation(data, output_file='reaction_animation.gif'):
    if not data:
        print("No data to animate.")
        return

    # Extract data arrays
    energies = [d['energy_ev'] for d in data]
    energies = np.array(energies)
    energies -= energies[0] # Relative to first frame
    
    # Smooth energy curve using simple moving average or spline if enough points
    # For short trajectories (e.g. 8 points), interpolation is better than smoothing window
    from scipy.interpolate import make_interp_spline
    
    num_points = len(energies)
    # Normalized Reaction Progress (0 to 1)
    reaction_progress = np.linspace(0, 1, num_points)
    
    # Create smooth curve for plotting line
    # Use more points for the smooth line
    reaction_progress_smooth = np.linspace(0, 1, 200)
    
    # Spline interpolation (k=3 for cubic, k=2 for quadratic)
    # Check if we have enough points for cubic spline (need > 3)
    k_order = 3 if num_points > 3 else (2 if num_points > 2 else 1)
    
    spl = make_interp_spline(reaction_progress, energies, k=k_order)
    energies_smooth = spl(reaction_progress_smooth)

    # Setup figure
    fig = plt.figure(figsize=(12, 6))
    
    # 3D Plot (Left)
    ax3d = fig.add_subplot(1, 2, 1, projection='3d')
    
    # Energy Plot (Right)
    ax_energy = fig.add_subplot(1, 2, 2)
    # Plot smooth curve for the background line
    ax_energy.plot(reaction_progress_smooth, energies_smooth, 'b-', label='Energy Profile', alpha=0.7)
    # Plot actual calculated points as faint markers
    ax_energy.plot(reaction_progress, energies, 'b.', alpha=0.3)
    
    ax_energy.set_xlabel(r'Reaction Progress $\xi$')
    ax_energy.set_ylabel('Energy (eV)')
    ax_energy.set_title('Reaction Energy Profile')
    energy_dot, = ax_energy.plot([], [], 'ro', markersize=10) # Moving dot

    # Pre-calculate bounds for 3D plot to keep camera steady
    all_coords = np.array([d['coordinates'] for d in data]) # Shape: (n_frames, n_atoms, 3)
    min_coord = all_coords.min(axis=(0, 1))
    max_coord = all_coords.max(axis=(0, 1))
    center = (min_coord + max_coord) / 2
    span = (max_coord - min_coord).max() / 2
    
    # Color map for charges
    cmap = plt.get_cmap('bwr') # Blue-White-Red (Negative-Neutral-Positive) or similar
    # Check charge range
    all_charges = []
    for d in data:
        all_charges.extend([a['charge'] for a in d['atoms']])
    max_charge = max(abs(min(all_charges)), abs(max(all_charges)))
    # Normalize charges to -max, +max or dynamic range
    # Use max_charge for symmetric range around 0 to highlight polarity
    # Or use min/max of actual data if not symmetric
    # Let's use symmetric range to keep 0 as white (neutral)
    norm = plt.Normalize(-max_charge, max_charge) 
    print(f"Charge color scale range: [-{max_charge:.3f}, {max_charge:.3f}]")

    def init():
        energy_dot.set_data([], [])
        return energy_dot,

    def update(frame_idx):
        frame_data = data[frame_idx]
        coords = np.array(frame_data['coordinates'])
        charges = [a['charge'] for a in frame_data['atoms']]
        symbols = frame_data['symbols']
        
        ax3d.clear()
        
        # Add colorbar (needs to be added once or managed carefully in animation)
        # Clearing ax3d removes it? No, colorbar is on figure usually. 
        # But if we add it every frame it stacks.
        # Let's add it outside update if possible, or use a mappable.
        
        # Plot atoms
        # Size based on element? (Optional)
        sizes = [100 if s != 'H' else 50 for s in symbols]
        
        # Color based on charge
        colors = cmap(norm(charges))
        
        sc = ax3d.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=sizes, c=colors, depthshade=True, edgecolors='k')
        
        # Draw bonds (simple distance criterion)
        # Naive bonding: < 1.6 Angstrom (C-C is ~1.5, C-H ~1.1, C-O ~1.4)
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < 1.6: 
                    ax3d.plot([coords[i,0], coords[j,0]], [coords[i,1], coords[j,1]], [coords[i,2], coords[j,2]], 'k-', linewidth=2, alpha=0.5)

        ax3d.set_xlim(center[0] - span, center[0] + span)
        ax3d.set_ylim(center[1] - span, center[1] + span)
        ax3d.set_zlim(center[2] - span, center[2] + span)
        ax3d.set_title(f'Frame {frame_idx} | Rel. Energy: {energies[frame_idx]:.2f} eV')
        
        # Update Energy Dot
        # Map frame index to reaction progress
        current_progress = frame_idx / (len(data) - 1) if len(data) > 1 else 0
        energy_dot.set_data([current_progress], [energies[frame_idx]]) # x and y must be sequences
        
        return ax3d, energy_dot

    # Add colorbar outside loop
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax3d, label='Atomic Charge (e)', shrink=0.5, aspect=10)

    # Slow down: interval 200ms -> 500ms (2 FPS)
    ani = FuncAnimation(fig, update, frames=len(data), init_func=init, blit=False, interval=500)
    
    print(f"Saving animation to {output_file}...")
    try:
        ani.save(output_file, writer='pillow', fps=2)
        print("Done!")
    except Exception as e:
        print(f"Error saving animation: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Animate reaction trajectory with DFT energy/charges.')
    parser.add_argument('--frames_dir', type=str, default="extracted_data_check", help='Directory containing XYZ frames')
    parser.add_argument('--last_frames', type=int, default=10, help='Number of last frames to process (default: 10)')
    parser.add_argument('--output', type=str, default='reaction_animation.gif', help='Output GIF filename')
    
    args = parser.parse_args()
    
    data = run_dft_on_frames(args.frames_dir, last_k_frames=args.last_frames)
    create_animation(data, output_file=args.output)
