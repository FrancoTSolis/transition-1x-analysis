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
python scripts/run_trajectory_analysis.py --h5_path Transition-State-Generation-Flow/dataset/transition1x/data/Transition1x.h5
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
    --h5_path Transition-State-Generation-Flow/dataset/transition1x/data/Transition1x.h5 \
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
