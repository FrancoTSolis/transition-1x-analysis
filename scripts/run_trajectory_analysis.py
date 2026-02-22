import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time
import random

# Add current directory to path to import extract_and_visualize
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import extract_and_visualize as ev

def run_atk_dft(xyz_file, script_path="scripts/run_dft.py"):
    """
    Runs the ATK DFT script for a given XYZ file.
    Returns the energy in Hartree.
    """
    cmd = ["atkpython", script_path, xyz_file]
    
    try:
        # Run command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse output for "DFT_ENERGY_HARTREE: <value>"
        for line in result.stdout.splitlines():
            if "DFT_ENERGY_HARTREE:" in line:
                energy_str = line.split(":")[1].strip()
                return float(energy_str)
        
        print(f"Warning: Could not find energy in output for {xyz_file}")
        print("Output:", result.stdout)
        return None
        
    except subprocess.CalledProcessError as e:
        print(f"Error running ATK for {xyz_file}:")
        print(e.stderr)
        return None
    except FileNotFoundError:
        print("Error: 'atkpython' not found. Make sure QuantumATK is installed and in your PATH.")
        return None

def process_trajectory(f, formula, rxn_id, args):
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

    # Run DFT for each frame
    calculated_energies = []
    frames_to_run = range(pos.shape[0])
    
    print(f"  Starting DFT calculations for {len(frames_to_run)} frames...")
    
    for i in frames_to_run:
        # Save XYZ
        xyz_file = os.path.join(xyz_dir, f"frame_{i:03d}.xyz")
        ev.save_xyz(xyz_file, z, pos[i], original_energies[i])
        
        if args.dry_run:
            # print(f"  Dry run: Skipping ATK for frame {i}")
            calculated_energies.append(original_energies[i]) 
            continue
            
        # Run ATK
        print(f"  Frame {i}/{len(frames_to_run)-1}...", end="", flush=True)
        t0 = time.time()
        energy = run_atk_dft(xyz_file)
        dt = time.time() - t0
        
        if energy is not None:
            print(f" Done ({dt:.2f}s). E={energy:.6f} Ha")
            calculated_energies.append(energy)
        else:
            print(" Failed.")
            calculated_energies.append(np.nan)

    # Plot Results
    plt.figure(figsize=(10, 6))
    
    orig_rel = original_energies - original_energies[0]
    
    valid_calc = np.array(calculated_energies)
    if not np.all(np.isnan(valid_calc)):
        first_valid_idx = np.where(~np.isnan(valid_calc))[0][0]
        calc_rel = valid_calc - valid_calc[first_valid_idx]
        calc_rel = calc_rel - calc_rel[first_valid_idx] + orig_rel[first_valid_idx]
    else:
        calc_rel = valid_calc

    plt.plot(range(len(original_energies)), orig_rel, 'o-', label='Original (wB97x/6-31G(d))', alpha=0.7)
    plt.plot(range(len(calculated_energies)), calc_rel, 's--', label='Calculated (ATK)', alpha=0.7)
    
    plt.xlabel('Frame Index')
    plt.ylabel('Relative Energy (Hartree)')
    plt.title(f'Energy Profile Comparison: {formula} - {rxn_id}')
    plt.legend()
    plt.grid(True)
    
    plot_path = os.path.join(reaction_dir, "energy_comparison.png")
    plt.savefig(plot_path)
    plt.close() # Close figure to free memory
    
    # Save raw data
    np.savetxt(os.path.join(reaction_dir, "calculated_energies.txt"), calculated_energies)
    print(f"  Analysis saved to {reaction_dir}")


def main():
    parser = argparse.ArgumentParser(description="Run trajectory analysis with QuantumATK")
    parser.add_argument("--h5_path", type=str, default="data/transition1x.h5", help="Path to h5 file")
    parser.add_argument("--output_dir", type=str, default="trajectory_analysis", help="Output directory")
    parser.add_argument("--formula", type=str, default="C2H2N2O", help="Chemical formula to select")
    parser.add_argument("--rxn_id", type=str, default="rxn2091", help="Reaction ID to select")
    parser.add_argument("--num_trajectories", type=int, default=1, help="Number of random trajectories to sample if formula/rxn_id not specified")
    parser.add_argument("--all", action="store_true", help="Run all trajectories in the dataset (overrides num_trajectories)")
    parser.add_argument("--dry_run", action="store_true", help="Don't actually run ATK, just extract and setup")
    
    args = parser.parse_args()
    
    try:
        f = ev.load_data(args.h5_path)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

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
            # Sample random trajectories
            num_to_sample = min(args.num_trajectories, total_reactions)
            trajectories_to_process = random.sample(all_reactions, num_to_sample)
            print(f"Randomly selected {num_to_sample} trajectories.")

    # Process each trajectory
    print(f"Processing {len(trajectories_to_process)} trajectories...")
    for i, (formula, rxn_id) in enumerate(trajectories_to_process):
        print(f"[{i+1}/{len(trajectories_to_process)}]", end=" ")
        process_trajectory(f, formula, rxn_id, args)

if __name__ == "__main__":
    main()
