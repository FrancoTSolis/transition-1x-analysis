# QuantumATK Trajectory Analysis Scripts

This directory contains scripts to extract transition state trajectories from the Transition1x dataset and perform Single Point Energy calculations using QuantumATK.

## 1. Setting up the Environment

To run these scripts, you need to set up a Python virtual environment that works with QuantumATK's python distribution (`atkpython`).

### Prerequisites
- QuantumATK X-2025.06 installed.
- Access to the dataset file `transition1x.h5`.

### Setup Instructions

1.  **Activate the QuantumATK Shell Environment**
    Run the following command in your terminal to set up the QuantumATK environment variables:
    ```bash
    source /media/francosolis/newdrive/qatk/tools/quantumatk/X-2025.06/atkpython/bin/activate_env
    ```
    Your prompt should change to include `(qatk)`.

2.  **Create the Virtual Environment**
    Create a new virtual environment at your desired location (`/media/francosolis/newdrive/atkpython_venv`) using the `atkpython` distribution. The `--system-site-packages` flag is **crucial** to allow access to QuantumATK modules.
    ```bash
    python -m venv --system-site-packages /media/francosolis/newdrive/atkpython_venv
    ```

3.  **Activate the Virtual Environment**
    ```bash
    source /media/francosolis/newdrive/atkpython_venv/bin/activate
    ```
    Your prompt should now look like `(atkpython_venv) (qatk) ...`.

4.  **Install Required Packages**
    Install the necessary Python libraries for data handling and plotting.
    ```bash
    pip install h5py matplotlib numpy
    ```
    *Note: `atkpython` comes with `numpy` and `scipy` pre-installed. Installing them again in the venv might upgrade them. If you encounter issues, try installing only missing packages like `h5py` (if not present) and `matplotlib`.*

## 2. Running the Analysis (Main Workflow)

The main script is `run_trajectory_analysis.py`. It handles data extraction, running QuantumATK DFT calculations, and plotting results.

### Basic Usage
Run analysis on a single random trajectory:
```bash
python scripts/run_trajectory_analysis.py --h5_path data/transition1x.h5
```

### Options

- **Sample Multiple Trajectories**:
  Run on 5 random trajectories:
  ```bash
  python scripts/run_trajectory_analysis.py --num_trajectories 5
  ```

- **Run All Trajectories**:
  (Warning: This might take a long time!)
  ```bash
  python scripts/run_trajectory_analysis.py --all
  ```

- **Run Specific Reaction**:
  Specify by formula and reaction ID:
  ```bash
  python scripts/run_trajectory_analysis.py --formula C2H2N2O --rxn_id rxn2091
  ```

- **Dry Run (Testing)**:
  Extract data and generate plots using original energies (skips actual DFT calculations):
  ```bash
  python scripts/run_trajectory_analysis.py --dry_run
  ```

### Output

Results are saved in the `trajectory_analysis/` directory (or whatever you specify with `--output_dir`).
Structure:
```
trajectory_analysis/
  └── Formula_RxnID/
      ├── xyz_frames/        # .xyz files for each frame
      ├── calculated_energies.txt
      └── energy_comparison.png
```

## 3. Running Individual Scripts

If you want to debug or run steps manually, you can use the individual helper scripts.

### Extracting and Visualizing Data (`extract_and_visualize.py`)

This script extracts a trajectory from the H5 file, saves XYZ frames, and plots the original energy profile without running any new calculations.

```bash
python scripts/extract_and_visualize.py \
    --h5_path data/transition1x.h5 \
    --formula C2H2N2O --rxn_id rxn2091 \
    --output_dir extracted_data_check
```

### Running DFT on a Single XYZ File (`run_dft.py`)

This is the QuantumATK script. It must be run with the `atkpython` executable (available after sourcing `activate_env`). It takes a single `.xyz` file and outputs the calculated energy.

```bash
# Ensure you are in the (qatk) environment
atkpython scripts/run_dft.py extracted_data_check/frame_000.xyz
```

*Note: This script prints the energy to stdout in the format `DFT_ENERGY_HARTREE: <value>`.*

## 4. Scripts Overview

- **`extract_and_visualize.py`**: Helper module for loading HDF5 data and visualizing frames. Can be run standalone to check data.
- **`run_dft.py`**: The actual script executed by `atkpython`. It performs the following:
  1.  **Reads XYZ**: Manually parses an XYZ file to create a `MoleculeConfiguration`. It maps element symbol strings (e.g., "O", "C") to QuantumATK `PeriodicTable` element objects.
  2.  **Calculator Setup**:
      -   **Basis Set**: Dynamically assigns `DoubleZetaPolarized` basis sets from the `GGABasis` module for each element present in the molecule (e.g., `GGABasis.Oxygen_DoubleZetaPolarized`). Falls back to other available basis sets if needed.
      -   **Exchange-Correlation**: Uses `GGA.PBE` as the functional. (Note: wB97X was originally requested but is not available in the standard namespace, so PBE is used as a robust standard).
  3.  **Calculation**: Runs a Self-Consistent Field (SCF) calculation to determine the electronic ground state.
  4.  **Energy and Charge Extraction**: 
      -   Calculates `TotalEnergy` and evaluates it to get the total energy in Hartree.
      -   Calculates `MullikenPopulation` to get net atomic charges.
  5.  **Output**: Prints a JSON object containing energy and per-atom charge data, framed by `JSON_OUTPUT_START` and `JSON_OUTPUT_END`. Also prints legacy `DFT_ENERGY_HARTREE` line.
- **`run_trajectory_analysis.py`**: The main driver script.
- **`animate_reaction.py`**: Creates a video animation of the reaction trajectory.
  -   **Inputs**: Directory of XYZ frames (e.g., `extracted_data_check`).
  -   **Process**:
      -   Subsamples frames (stride ~20-30 frames total) for performance.
      -   Runs `run_dft.py` on each selected frame to get energy and charges.
      -   Caches results in `trajectory_data.json` to avoid re-calculation.
  -   **Output**: Generates `reaction_animation.gif` with:
      -   **Left Panel**: 3D visualization of the molecule, atoms colored by partial charge (Blue=Positive, Red=Negative).
      -   **Right Panel**: Energy profile with a moving marker indicating the current frame.
  -   **Usage**: `python scripts/animate_reaction.py --last_frames 8`

## 5. ANI-2x Scripts (No QuantumATK Required)

A parallel set of scripts uses the **ANI-2x neural network potential** (`torchani`) instead of QuantumATK DFT. These run with standard Python and are significantly faster since energies are computed in batch on GPU/CPU without spawning subprocesses.

### Prerequisites (ANI)

```bash
pip install torch torchani h5py matplotlib numpy scipy
```

ANI-2x supports a limited set of elements: **H, C, N, O, F, S, Cl**.

### Running ANI Energy on a Single XYZ File (`run_ani.py`)

The ANI counterpart of `run_dft.py`. No `atkpython` needed.

```bash
python scripts/run_ani.py extracted_data_check/frame_000.xyz
# With GPU:
python scripts/run_ani.py extracted_data_check/frame_000.xyz --device cuda
```

Output format is compatible with `run_dft.py` (JSON block + `ANI_ENERGY_HARTREE: <value>`).

### Running ANI Trajectory Analysis (`run_trajectory_analysis_ani.py`)

The ANI counterpart of `run_trajectory_analysis.py`. Computes all frame energies in a single batch.

```bash
# Specific reaction
python scripts/run_trajectory_analysis_ani.py --formula C2H2N2O --rxn_id rxn2091

# Dry run
python scripts/run_trajectory_analysis_ani.py --dry_run

# Use GPU
python scripts/run_trajectory_analysis_ani.py --device cuda
```

Output is saved to `trajectory_analysis_ani/` by default, with the same directory structure as the DFT version plus an additional `summary_ani.json`.

### Animating with ANI Energies (`animate_reaction_ani.py`)

The ANI counterpart of `animate_reaction.py`. Atoms are colored by element (CPK convention) since ANI does not provide partial charges.

```bash
python scripts/animate_reaction_ani.py \
    --frames_dir trajectory_analysis_ani/C2H2N2O_rxn2091/xyz_frames \
    --last_frames 8
```

## 6. Ground-Truth Animation (`animate_reaction_gt.py`)

Animates the trajectory using the **ground-truth wB97x/6-31G(d) energies** that are already embedded in each XYZ file's comment line. No DFT or ANI computation is performed — only `numpy`, `matplotlib`, and `scipy` are required.

```bash
python scripts/animate_reaction_gt.py \
    --frames_dir trajectory_analysis_ani/C2H2N2O_rxn2091/xyz_frames \
    --last_frames 8
```

Output: `reaction_animation_gt.gif`.

## 7. The `--include_reactant_and_product` Flag

All three animation scripts (`animate_reaction.py`, `animate_reaction_ani.py`, `animate_reaction_gt.py`) support this flag. It is **on by default**.

When enabled, the selected frames follow this pattern (for `--last_frames K`):

```
[0,  -K, -(K-1), ..., -1,  K+1]
 ^   \_________________/    ^
 |     last K frames        |
 reactant                 product
```

- **Frame 0**: the reactant geometry (first frame of the trajectory).
- **Frames -K to -1**: the last K frames of the trajectory.
- **Frame K+1**: the product geometry.

To disable and only use the last K frames (original behaviour):

```bash
python scripts/animate_reaction_gt.py --last_frames 8 --no-include_reactant_and_product
```

## 8. Energy Units — Important Notes

The three methods output energies in **different units**. Be aware of this when comparing:

| Source | Method | Native Unit | Typical Value (C2H2N2O) |
|---|---|---|---|
| **Ground truth** | wB97x/6-31G(d) from Transition1x H5 | **eV** | -7130 eV |
| **ANI-2x** | `torchani.models.ANI2x` | **Hartree** | -262 Ha |
| **DFT (ATK)** | QuantumATK `LCAOCalculator` (GGA.PBE) | **Hartree** | ~-262 Ha |

Conversion: **1 Hartree = 27.211386 eV**.

The ground-truth energies stored in the Transition1x dataset (`wB97x_6-31G(d).energy` field in the H5 file) are in **eV**. These are written into the XYZ comment line as `Energy=<value>` when frames are extracted by `extract_and_visualize.py`. The `animate_reaction_gt.py` script reads these values directly — no conversion is needed.

ANI-2x and QuantumATK both output in **Hartree**. The animation scripts (`animate_reaction_ani.py`, `animate_reaction.py`) convert to eV internally before plotting.

*Note: The key `original_energies_hartree` in `summary_ani.json` is a misnomer — those values are actually in eV (from the Transition1x dataset).*

## 9. 3-Way Energy Comparison (`compare_energy_methods.py`)

Plots DFT, ground truth, and ANI-2x energy curves side by side for the same selected frames. Ground truth is in the middle panel. The energy barrier (ΔE‡ = E_TS − E_reactant) is annotated on all three panels, where the TS frame is the one with the highest ground-truth energy among the selected frames. On the DFT and ANI panels, the ground-truth curve is overlaid as a dashed reference and the per-frame error is shaded with red hatched lines.

Requires: numpy, matplotlib (no torch or atkpython).

### Usage

```bash
python scripts/compare_energy_methods.py \
    --xyz_dir trajectory_analysis_ani/C2H2N2O_rxn2091/xyz_frames \
    --ani_energies trajectory_analysis_ani/C2H2N2O_rxn2091/calculated_energies_ani.txt \
    --dft_cache trajectory_data.json \
    --last_frames 8
```

### Options

- `--xyz_dir` (required): directory with `frame_*.xyz` files containing ground-truth energies in the comment line.
- `--ani_energies`: path to `calculated_energies_ani.txt` (Hartree, one value per line).
- `--dft_cache`: path to `trajectory_data.json` (JSON cache from `animate_reaction.py`, energies in eV).
- `--dft_energies`: alternative to `--dft_cache` — plain-text file of DFT energies in Hartree.
- `--last_frames N`: number of last frames to include (default: 8).
- `--include_reactant_and_product` / `--no-include_reactant_and_product`: include frame 0 (reactant) and frame N+1 (product) in the selection (default: on).
- `--output`: output image filename (default: `energy_comparison_3way.png`).

### Output

`energy_comparison_3way.png` with three panels:

| Left | Middle | Right |
|---|---|---|
| DFT — GGA.PBE | Ground Truth — ωB97x/6-31G(d) | ANI-2x |
| + red hatched error vs GT | energy barrier annotated | + red hatched error vs GT |

## 10. Full Scripts Overview

| Script | Requires | Purpose |
|---|---|---|
| `extract_and_visualize.py` | numpy, h5py, matplotlib | Extract trajectories from H5, save XYZ frames, plot energy profiles |
| `run_dft.py` | **atkpython** | Single-point DFT energy via QuantumATK (GGA.PBE) |
| `run_ani.py` | torch, torchani | Single-point energy via ANI-2x neural network potential |
| `run_trajectory_analysis.py` | **atkpython**, h5py, numpy, matplotlib | Batch trajectory analysis with QuantumATK DFT |
| `run_trajectory_analysis_ani.py` | torch, torchani, h5py, numpy, matplotlib | Batch trajectory analysis with ANI-2x |
| `animate_reaction.py` | **atkpython**, numpy, matplotlib, scipy | Animate trajectory with DFT energies + Mulliken charges |
| `animate_reaction_ani.py` | torch, torchani, numpy, matplotlib, scipy | Animate trajectory with ANI-2x energies |
| `animate_reaction_gt.py` | numpy, matplotlib, scipy | Animate trajectory with ground-truth wB97x/6-31G(d) energies (no model needed) |
| `compare_energy_methods.py` | numpy, matplotlib | 3-way side-by-side comparison of DFT / GT / ANI energy curves |
