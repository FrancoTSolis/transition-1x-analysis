# qchem_pipeline

A CLI tool that automates the full Q-Chem workflow for transition-state
benchmark calculations: **generate** input files from XYZ coordinates, **run**
them locally or on a SLURM cluster, and **parse** the outputs into structured
energy/amplitude summaries.

```
XYZ files ──→ generate ──→ .in files ──→ run ──→ .out files ──→ parse ──→ JSON/CSV
                                         │
                                    local or SLURM
```

---

## Setup

```bash
# Only dependency is PyYAML
pip install pyyaml

# Or from the requirements file
pip install -r qchem_pipeline/requirements.txt
```

The tool is invoked as a Python module from the `ciqc_collab/` directory:

```bash
cd ciqc_collab/
python -m qchem_pipeline --help
```

> Throughout this document, `python` refers to your environment's Python.
> On this machine: `/xuanwu-tank/east/fts/conda_envs/t1x/bin/python`

---

## Quick Start

### One-command pipeline (local, no SLURM)

```bash
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/ja025633n_sto3g/ \
    --basis STO-3G --method "ccsd(t)" --mode local
```

This creates `runs/ja025633n_sto3g/` with three subdirectories:
- `inputs/` — generated `.in` files
- `outputs/` — Q-Chem `.out` and `.fchk` files
- `results.json` and `results.csv` — parsed energy summary

### One-command pipeline (SLURM)

```bash
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/ja025633n_631gd/ \
    --basis "6-31G*" --method "ccsd(t)" --mem 24000 \
    --mode slurm --account MY_ALLOCATION --time 08:00:00
```

After jobs complete, parse separately:

```bash
python -m qchem_pipeline parse runs/ja025633n_631gd/outputs/ \
    --json runs/ja025633n_631gd/results.json \
    --csv  runs/ja025633n_631gd/results.csv
```

### Dry run (preview without executing)

Append `--dry-run` to any `run` or `pipeline` command to see what would happen:

```bash
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/test/ \
    --basis STO-3G --mode local --dry-run
```

---

## Commands

### `generate` — Create Q-Chem input files from XYZ

```
python -m qchem_pipeline generate <xyz_dir> <output_dir> [options]
```

| Argument | Description |
|----------|-------------|
| `xyz_dir` | Directory containing `.xyz` files |
| `output_dir` | Where to write generated `.in` files |
| `--method`, `-m` | Method override: `auto`, `ccsd`, `ccsd(t)`, or `mp2` |
| `--basis`, `-b` | Basis set override: `STO-3G`, `6-31G*`, `cc-pVTZ`, etc. |
| `--mem` | `MEM_TOTAL` in MB (omit to let Q-Chem decide) |

**Behavior:**
- Reads each `.xyz` file, extracts coordinates (skipping the atom-count and
  comment lines).
- Looks up charge/multiplicity from the molecule registry (see Configuration).
- **`auto` method** (default): molecules with ≤40 atoms get `ccsd(t)`, larger
  ones get `mp2`. This is the recommended setting for mixed-size datasets.
- If a specific method like `ccsd(t)` is requested but the molecule has >40
  atoms, it is still downgraded to MP2 with a warning.
- Generates one `.in` file per molecule, named
  `{molecule}_{method}_{basis}.in`.

**Examples:**

```bash
# Small basis for pipeline validation
python -m qchem_pipeline generate \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/sto3g/inputs/ \
    --basis STO-3G --method "ccsd(t)"

# Production basis with memory allocation
python -m qchem_pipeline generate \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/631gd/inputs/ \
    --basis "6-31G*" --method "ccsd(t)" --mem 24000

# Force MP2 for all molecules
python -m qchem_pipeline generate \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/mp2_test/inputs/ \
    --method mp2 --basis STO-3G
```

---

### `run` — Execute Q-Chem jobs

```
python -m qchem_pipeline run <input_dir> [options]
```

| Argument | Description |
|----------|-------------|
| `input_dir` | Directory containing `.in` files |
| `--mode` | `local` (default), `slurm`, or `slurm-generate` |
| `--output-dir`, `-o` | Where to write `.out` files (default: same as input) |
| `--nprocs`, `-np` | Number of CPU threads for Q-Chem |
| `--qchem-cmd` | Path to `qchem` executable (default: `qchem`) |
| `--account` | SLURM account/allocation (required for `slurm` mode) |
| `--time` | SLURM wall time, e.g. `08:00:00` |
| `--dry-run` | Print commands without executing |

**Three modes:**

#### `--mode local` (default)

Runs each `.in` file sequentially via `qchem` subprocess:

```bash
python -m qchem_pipeline run runs/sto3g/inputs/ --mode local --nprocs 4
```

#### `--mode slurm`

Generates a `.slm` SLURM batch script for each `.in` file and submits it
via `sbatch`:

```bash
python -m qchem_pipeline run runs/631gd/inputs/ \
    --mode slurm \
    --account MY_ALLOCATION \
    --time 04:00:00 \
    --nprocs 32
```

#### `--mode slurm-generate`

Generates `.slm` scripts without submitting — useful when you want to inspect
or edit them before running `sbatch` manually:

```bash
python -m qchem_pipeline run runs/631gd/inputs/ \
    --mode slurm-generate \
    --account MY_ALLOCATION \
    --time 08:00:00

# Inspect, then submit manually
cat runs/631gd/inputs/molecule_name.slm
sbatch runs/631gd/inputs/molecule_name.slm
```

---

### `parse` — Extract results from Q-Chem output

```
python -m qchem_pipeline parse <output_dir> [options]
```

| Argument | Description |
|----------|-------------|
| `output_dir` | Directory containing `.out` files |
| `--json` | Export full results to a JSON file |
| `--csv` | Export energy summary to a CSV file |

**What it extracts from each `.out` file:**

| Field | Description |
|-------|-------------|
| SCF energy | Hartree-Fock energy |
| MP2 energy | MP2 total energy |
| CCSD energy | CCSD total energy |
| CCSD(T) energy | CCSD(T) total energy (if triples were computed) |
| CCSD correlation energy | The correlation correction added to HF |
| T1², T2² | Coupled-cluster amplitude diagnostics |
| Leading amplitudes | Largest t₁ and t₂ excitation coefficients |
| SCF/CC convergence | Whether each stage converged |
| .fchk existence | Whether the formatted checkpoint file was generated |
| Wall time | Total job time |

**Example output:**

```
Molecule                                            OK   Method      Basis            E(HF)           E(MP2)          E(CCSD)       E(CCSD(T))     T1²       Time
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
F2_aug_dz                                            ✓  ccsd(T) aug-cc-pVdz    -198.69867529    -199.12680924    -199.13419803    -199.14706241  0.0055 13.34s(wall)

  1/1 succeeded
```

**Examples:**

```bash
# Print summary to terminal
python -m qchem_pipeline parse runs/sto3g/outputs/

# Export to both formats
python -m qchem_pipeline parse runs/sto3g/outputs/ \
    --json results.json --csv results.csv
```

---

### `pipeline` — Full end-to-end workflow

```
python -m qchem_pipeline pipeline <xyz_dir> <work_dir> [options]
```

Combines `generate` → `run` → `parse` into a single command.

| Argument | Description |
|----------|-------------|
| `xyz_dir` | Directory containing `.xyz` files |
| `work_dir` | Working directory (will contain `inputs/`, `outputs/`, results) |
| `--mode` | `local` (default) or `slurm` |
| `--method`, `-m` | Method override |
| `--basis`, `-b` | Basis set override |
| `--mem` | `MEM_TOTAL` in MB |
| `--nprocs`, `-np` | CPU threads |
| `--qchem-cmd` | Q-Chem executable path |
| `--account` | SLURM account (required for `--mode slurm`) |
| `--time` | SLURM wall time |
| `--dry-run` | Preview without executing |

**Directory structure created:**

```
work_dir/
├── inputs/          # Generated .in files
├── outputs/         # Q-Chem .out and .fchk files
├── results.json     # Full parsed results
└── results.csv      # Flat energy summary
```

> In SLURM mode, the pipeline stops after submitting jobs. Run `parse`
> separately once jobs complete.

---

### `status` — Check job progress

```
python -m qchem_pipeline status <output_dir>
```

- If `.out` files exist: parses them and prints the summary table.
- If no outputs yet and SLURM is available: queries `squeue` for your jobs.
- Otherwise: reports how many input files are waiting.

```bash
python -m qchem_pipeline status runs/631gd/outputs/
```

---

## Configuration

The pipeline looks for `pipeline_config.yaml` in the current directory or any
parent directory. CLI flags override config values.

```yaml
qchem:
  method: "auto"                # auto, ccsd, ccsd(t), or mp2
  basis: "STO-3G"              # Basis set
  mem_total: null               # MEM_TOTAL in MB (null = auto)
  unrestricted: true            # UHF for transition states
  internal_stability: true      # Check SCF orbital stability
  scf_algorithm: "DIIS_GDM"    # DIIS_GDM or GDM
  max_scf_cycles: 200
  thresh: 15                    # Integral threshold
  scf_convergence: 10          # SCF convergence threshold
  cc_convergence: 7             # CC convergence threshold
  n_frozen_core: "FC"           # FC = standard frozen core
  gui: 2                        # Generate .fchk for MO coefficients

slurm:
  account: null                 # Your SLURM allocation
  time: "08:00:00"
  nodes: 1
  ntasks: 1
  cpus_per_task: 18
  constraint: "cpu"
  qos: "regular"
  modules: [qchem]

# Per-molecule charge/multiplicity overrides (matched by filename prefix)
charge_mult_overrides:
  "02_Cationic_Aza_Cope": [1, 1]
```

**Defaults without a config file:** charge=0, multiplicity=1, STO-3G,
auto method (ccsd(t) for ≤40 atoms, mp2 for larger), unrestricted HF,
DIIS_GDM, GUI=2.

---

## Recipes

### Running the Epoxide / Sulfur Ylide dataset (ja025633n_s1)

12 molecules, 16–36 atoms each. All neutral singlets (charge=0, mult=1).
All are small enough for CCSD(T) with STO-3G; the 36-atom molecules may
need MP2 with larger bases.

```bash
# --- STO-3G validation (all CCSD(T), finishes in minutes) ---
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/ja025633n_sto3g/ \
    --basis STO-3G --method "ccsd(t)" --mode local

# --- Production basis, local ---
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/ja025633n_631gd/ \
    --basis "6-31G*" --mem 24000 --mode local --nprocs 8

# --- Production basis, SLURM ---
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/ja025633n_ccpvtz/ \
    --basis "cc-pVTZ" --mem 100000 \
    --mode slurm --account MY_ACCT --time 12:00:00

# Then parse after jobs finish:
python -m qchem_pipeline parse runs/ja025633n_ccpvtz/outputs/ \
    --json runs/ja025633n_ccpvtz/results.json \
    --csv  runs/ja025633n_ccpvtz/results.csv
```

### Running the Sigmatropic Rearrangement dataset (ja4c17174_si_001)

6 molecules with a wide size range: 14–93 atoms. Molecules 01–02 are small
(14–16 atoms, CCSD(T) feasible), molecules 03–06 are large (82–93 atoms,
MP2 only). Molecule 02 is a **cation** (charge=+1), handled automatically
via `pipeline_config.yaml`.

```bash
# --- Auto method (recommended): CCSD(T) for small, MP2 for large ---
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/ja4c17174_sto3g/ \
    --basis STO-3G --mode local

# This is equivalent (auto is the default):
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/ja4c17174_sto3g/ \
    --method auto --basis STO-3G --mode local

# Output will show:
#   01_..._ccsdt_STO-3G.in  (14 atoms, ccsd(t), charge=0, mult=1)
#   02_..._ccsdt_STO-3G.in  (16 atoms, ccsd(t), charge=1, mult=1)
#   03_..._mp2_STO-3G.in    (93 atoms, mp2,     charge=0, mult=1)
#   04_..._mp2_STO-3G.in    (93 atoms, mp2,     charge=0, mult=1)
#   05_..._mp2_STO-3G.in    (82 atoms, mp2,     charge=0, mult=1)
#   06_..._mp2_STO-3G.in    (82 atoms, mp2,     charge=0, mult=1)

# --- Force MP2 for all (e.g. for a uniform comparison) ---
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/ja4c17174_mp2_sto3g/ \
    --method mp2 --basis STO-3G --mode local

# --- Production basis via SLURM ---
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/ja4c17174_631gd/ \
    --basis "6-31G*" --mem 48000 \
    --mode slurm --account MY_ACCT --time 08:00:00
```

### Step-by-step SLURM workflow (either dataset)

For maximum control — generate, inspect, submit, parse separately:

```bash
# Step 1: Generate inputs
python -m qchem_pipeline generate \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/production/inputs/ \
    --basis "cc-pVTZ" --method auto --mem 100000

# Step 2: Generate SLURM scripts (inspect before submitting)
python -m qchem_pipeline run runs/production/inputs/ \
    --mode slurm-generate \
    -o runs/production/outputs/ \
    --account MY_ACCT --time 12:00:00 --nprocs 32

# Step 3: Submit all
for f in runs/production/outputs/*.slm; do sbatch "$f"; done

# Step 4: Monitor
python -m qchem_pipeline status runs/production/outputs/

# Step 5: Parse when done
python -m qchem_pipeline parse runs/production/outputs/ \
    --json runs/production/results.json \
    --csv  runs/production/results.csv
```

### Other recipes

```bash
# Dry-run: preview what would happen without executing
python -m qchem_pipeline pipeline \
    coordinates_xyz_output/ja4c17174_si_001/ \
    runs/test/ \
    --basis STO-3G --mode local --dry-run

# Custom config file
python -m qchem_pipeline -c my_custom_config.yaml pipeline \
    coordinates_xyz_output/ja025633n_s1/ \
    runs/custom/ \
    --mode local
```

---

## Molecule Charge/Multiplicity

Most molecules default to **charge=0, multiplicity=1** (neutral singlet).

Exceptions are specified in `pipeline_config.yaml` under
`charge_mult_overrides`. The key is matched as a **filename prefix**:

```yaml
charge_mult_overrides:
  "02_Cationic_Aza_Cope": [1, 1]    # cation, singlet
  "some_radical":         [0, 2]    # neutral, doublet
```

If your XYZ file is named `02_Cationic_Aza_Cope_Rearrangement_Transition_State.xyz`,
the prefix `02_Cationic_Aza_Cope` matches and sets charge=1, mult=1.

---

## Output Formats

### JSON (`--json`)

Full structured output including per-molecule SCF/CC convergence details,
leading amplitudes, and metadata:

```json
[
  {
    "filename": "molecule.out",
    "success": true,
    "method": "ccsd(t)",
    "basis": "STO-3G",
    "scf": { "converged": true, "energy": -198.698, "cycles": 10 },
    "mp2": { "energy": -199.126 },
    "cc": {
      "converged": true,
      "ccsd_energy": -199.134,
      "ccsd_t_energy": -199.147,
      "t1_squared": 0.0055,
      "t2_squared": 0.1179
    },
    "leading_amplitudes": [ ... ]
  }
]
```

### CSV (`--csv`)

Flat table for quick spreadsheet analysis:

```
molecule,success,method,basis,natoms,scf_energy,mp2_energy,ccsd_energy,ccsd_t_energy,...
```

---

## Python API

The modules can also be used directly in Python scripts:

```python
from qchem_pipeline.config import ProjectConfig
from qchem_pipeline.generate import generate_inputs, read_xyz
from qchem_pipeline.run import run_local, run_local_batch
from qchem_pipeline.parse import parse_output, parse_batch, print_summary

# Load config
config = ProjectConfig.from_yaml("pipeline_config.yaml")

# Generate inputs
generated = generate_inputs(
    "coordinates_xyz_output/ja025633n_s1/",
    "runs/test/inputs/",
    config,
    basis_override="STO-3G",
)

# Parse a single output
result = parse_output("runs/test/outputs/molecule.out")
print(result.scf.energy)
print(result.cc.ccsd_t_energy)
print(result.cc.t1_squared)

# Parse a whole directory
results = parse_batch("runs/test/outputs/")
print_summary(results)
```
