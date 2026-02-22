# Q-Chem CCSD/MP2 Pipeline Guide for LUCJ Ansatz Initialization

> **Audience:** AI/quantum-computing researcher who needs high-precision classical
> amplitudes (t₁, t₂) and MO coefficients from CCSD(T) and MP2 calculations on
> transition-state geometries, to seed a LUCJ (Locally-Unitary Coupled-Cluster
> Jastrow) ansatz on a quantum circuit.

---

## Table of Contents

1. [The Big Picture: Why We're Doing This](#1-the-big-picture)
2. [Background: What Are CCSD, CCSD(T), and MP2?](#2-background)
3. [Your Molecule Inventory](#3-molecule-inventory)
4. [Anatomy of the Q-Chem Templates](#4-templates)
5. [End-to-End Pipeline](#5-pipeline)
6. [Step-by-Step: Generating Input Files from XYZ](#6-generating-inputs)
7. [Understanding the Output](#7-understanding-output)
8. [Extracting Amplitudes and MO Coefficients](#8-extracting-data)
9. [Practical Tips and Troubleshooting](#9-tips)
10. [Recommended Execution Plan](#10-execution-plan)
11. [Questions to Clarify with Collaborator](#11-questions)

---

## 1. The Big Picture

```
Literature TS geometries (XYZ)
        │
        ▼
┌──────────────────────┐
│  Q-Chem: HF → CCSD  │  ← classical computation
│  (or HF → MP2)       │
└──────────┬───────────┘
           │
    ┌──────┴──────┐
    │             │
    ▼             ▼
 Energies    Amplitudes + MO coefficients
 (benchmarks)  (t₁, t₂ from CC; MO coeffs from HF)
                  │
                  ▼
        ┌─────────────────────┐
        │ Quantum Circuit:    │
        │ LUCJ Ansatz Init    │
        │ (warm-start)        │
        └─────────────────────┘
```

**Goal:** We are *not* optimizing geometries — the geometries are already optimized
(taken from literature). We run **single-point energy** calculations to get:

1. **Benchmark energies** at CCSD(T) and MP2 levels of theory.
2. **t₁ and t₂ amplitudes** (the coupled-cluster excitation coefficients) — these
   are used to initialize the unitary coupled-cluster parameters on the quantum
   circuit.
3. **MO (molecular orbital) coefficients** — the transformation matrix from
   atomic orbitals (AOs) to molecular orbitals (MOs), needed to define the
   second-quantized Hamiltonian in the MO basis.

---

## 2. Background: What Are CCSD, CCSD(T), and MP2?

If you're coming from an AI/CS background, think of these as progressively more
expensive (but more accurate) ways to solve the electronic Schrödinger equation.

### The Hierarchy (Jacob's Ladder of Electron Correlation)

| Method | Accuracy | Scaling | Analogy |
|--------|----------|---------|---------|
| **HF** (Hartree-Fock) | Baseline — mean-field, no correlation | O(N⁴) | Like a single linear layer |
| **MP2** (2nd-order Møller-Plesset) | Captures ~80-90% of correlation energy | O(N⁵) | A perturbative correction — like a first-order Taylor expansion around HF |
| **CCSD** (Coupled Cluster Singles & Doubles) | Very accurate for single-ref systems | O(N⁶) | An exponential ansatz that systematically includes single and double excitations |
| **CCSD(T)** | "Gold standard" of quantum chemistry | O(N⁷) | CCSD + a perturbative triples correction |

### What Actually Happens in a Q-Chem Run

Every calculation follows the same two-stage flow:

```
Stage 1: Self-Consistent Field (SCF) = Hartree-Fock
  ┌──────────────────────────────────────────┐
  │ Iteratively solve for the MOs (orbitals) │
  │ that minimize the mean-field energy.     │
  │                                          │
  │ Output: HF energy, MO coefficients,      │
  │         orbital energies                  │
  └────────────────────┬─────────────────────┘
                       │
Stage 2: Post-HF correlation (CCSD or MP2)
  ┌──────────────────────────────────────────┐
  │ Use the HF orbitals as a basis and       │
  │ compute electron correlation on top.     │
  │                                          │
  │ Output: Correlation energy, t₁/t₂ amps, │
  │         total energy                     │
  └──────────────────────────────────────────┘
```

### Key Concepts for the AI-Minded

- **Basis set**: The set of mathematical functions used to expand the MOs. Larger
  basis sets = more accurate but exponentially more expensive. Think of it as the
  "resolution" of the calculation.
  - `STO-3G`: minimal basis (~3 functions per atom). Fast, inaccurate. **Use for
    testing your pipeline.**
  - `6-31G(d)`: small-medium. Reasonable for qualitative results.
  - `cc-pVDZ`, `cc-pVTZ`, `def2-TZVPD`: larger, production-quality bases.

- **Frozen core**: Inner-shell electrons (like 1s on carbon) are chemically inert.
  "Freezing" them (excluding from correlation) saves cost with negligible accuracy
  loss. Q-Chem does this by default (`N_FROZEN_CORE = FC`).

- **Charge & multiplicity**: Every molecule has a net charge (0 for neutral, +1
  for cation, etc.) and a spin multiplicity (1 = singlet = all electrons paired,
  2 = doublet = one unpaired electron, etc.). These must be specified correctly.

- **Unrestricted vs. Restricted HF**: For singlet states, restricted HF (RHF)
  forces alpha and beta electrons into the same spatial orbitals. Unrestricted
  (UHF) lets them differ, which is important for transition states that may have
  partial bond-breaking character where the restricted solution is unstable.

- **t₁ and t₂ amplitudes**: In coupled-cluster theory, the wavefunction is
  written as |Ψ⟩ = e^T |HF⟩, where T = T₁ + T₂ + .... T₁ encodes single
  excitations (one electron jumps from occupied → virtual orbital), T₂ encodes
  double excitations (two electrons jump simultaneously). The numerical
  coefficients of these excitations are the **t amplitudes** — these are exactly
  what you need to initialize the LUCJ ansatz.

---

## 3. Your Molecule Inventory

### Dataset 1: Sigmatropic Rearrangements (ja4c17174)

Source: M06-2X/6-311+G(d,p) optimized geometries, original ref uses CCSD(T)/cc-pVTZ single-point energies.

| # | Name | Atoms | Elements | Charge | Mult | Tier |
|---|------|-------|----------|--------|------|------|
| 01 | Claisen Rearrangement TS | 14 | C,H,O | 0 | 1 | CCSD(T) feasible |
| 02 | Cationic Aza-Cope TS | 16 | C,H,N | **1** | 1 | CCSD(T) feasible |
| 03 | TS-1 ketenimine exo-lone pair | **93** | C,H,N,O,Si | 0 | 1 | MP2 only |
| 04 | TS-SI-1 ketenimine endo-lone pair | **93** | C,H,N,O,Si | 0 | 1 | MP2 only |
| 05 | TS-2 formaldimine 5 exo-lone pair | **82** | C,H,N,O,Si | 0 | 1 | MP2 only |
| 06 | TS-SI-2 formaldimine 5 endo-lone pair | **82** | C,H,N,O,Si | 0 | 1 | MP2 only |

> **Important**: Molecules 03-06 are very large (82-93 atoms). CCSD(T) is not
> feasible for these. Use MP2, or discuss with the collaborator about active-space
> methods. Even MP2 with a large basis may be expensive.

> **Note on charge/multiplicity**: Molecule 02 is explicitly cationic (charge=+1)
> per the comment line. The others are neutral singlets. **Verify charge/mult for
> molecules 03-06 with the collaborator** — they contain Si and could have unusual
> electronic structure.

### Dataset 2: Epoxide Formation from Sulfur Ylides (ja025633n)

Source: B3LYP/6-31+G* (or 6-31G*) optimized geometries. Reference MP2/6-311+G** energies available in comment lines.

| # | Name | Atoms | Elements | Charge | Mult | Tier |
|---|------|-------|----------|--------|------|------|
| 01 | Torsional Rotation TS (7) | 16 | C,H,O,S | 0 | 1 | CCSD(T) feasible |
| 02 | Transoid Elimination TS (9) | 16 | C,H,O,S | 0 | 1 | CCSD(T) feasible |
| 03 | Cis Elimination TS (12) | 16 | C,H,O,S | 0 | 1 | CCSD(T) feasible |
| 04 | Syn Cisoid Addition TS (13) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 05 | Anti Cisoid Addition TS (14) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 06 | Anti Cisoid Addition TS rotamer (15) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 07 | Syn Transoid Addition TS | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 08 | Syn Torsional Rotation TS (18) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 09 | Anti Torsional Rotation TS (19) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 10 | Syn SN2 Elimination TS (25) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 11 | Anti SN2 Elimination TS (26) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |
| 12 | Cis Elimination TS (24) | 36 | C,H,O,S | 0 | 1 | CCSD(T) or MP2 |

> The 16-atom molecules are small enough for CCSD(T) with any basis. The 36-atom
> molecules are borderline — CCSD(T)/STO-3G should work; CCSD(T) with a large
> basis may require significant resources or MP2 as a fallback.

---

## 4. Anatomy of the Q-Chem Templates

### 4.1 CCSD Template (`ccsd_template.in`)

```
$molecule
TODO_CHARGE TODO_MULT
TODO_MOL
$end

$rem
    METHOD          ccsd
    BASIS           TODO_BASIS
    MEM_TOTAL       mem_opt
    UNRESTRICTED    hf_opt
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   scf_opt
    MAX_SCF_CYCLES  scf_cycle_opt
    THRESH          thresh_opt
    SCF_CONVERGENCE conv_opt
    CC_CONVERGENCE  cc_conv_opt
$end
```

**Line-by-line explanation:**

| Variable | What it does | Recommended value |
|----------|-------------|-------------------|
| `TODO_CHARGE` | Net charge of the molecule | 0 for most molecules; 1 for Aza-Cope |
| `TODO_MULT` | Spin multiplicity (2S+1) | 1 for closed-shell singlets |
| `TODO_MOL` | Cartesian coordinates (paste from XYZ, **without** the atom-count and comment lines) | Copy lines 3+ from the .xyz file |
| `METHOD ccsd` | Request CCSD calculation. Q-Chem will also compute CCSD(T) automatically when you add `CCSD_T` or use `METHOD ccsd(t)` | Keep as-is, or change to `ccsd(t)` for the triples correction |
| `TODO_BASIS` | Basis set name | Start with `STO-3G`, then `6-31G*`, then target basis |
| `mem_opt` | Total memory in MB (1000 = 1 GB) | Omit for STO-3G; `24000`-`300000` for production runs |
| `UNRESTRICTED` | Use UHF instead of RHF | `TRUE` for transition states (recommended); `FALSE` for "easy" systems (halves memory) |
| `INTERNAL_STABILITY` | Check if the SCF solution is a true minimum (not a saddle point) in orbital space | Always `TRUE` — critical for transition states |
| `SCF_ALGORITHM` | Algorithm for SCF convergence | `DIIS_GDM` (try first); `GDM` (if DIIS_GDM fails) |
| `MAX_SCF_CYCLES` | Max iterations for HF convergence | `200` is generous; `100` is default |
| `THRESH` | Integral screening threshold (10^-THRESH) | `15` (tighter than default; ensures precision) |
| `SCF_CONVERGENCE` | SCF converged when DIIS error < 10^-conv | `10` (very tight; lower to 8 if near-convergence stalls) |
| `CC_CONVERGENCE` | CC converged when amplitude change < 10^-conv | `7` (good default) |

> **Note:** The template says `METHOD ccsd`, not `ccsd(t)`. To get the perturbative
> triples correction (T), change it to `METHOD ccsd(t)`. Q-Chem will then report
> CCSD, (T) correction, and CCSD(T) total energy separately.

### 4.2 MP2 Template (`mp2_template.in`)

```
$molecule
TODO_CHARGE TODO_MULT
TODO_MOL
$end

$rem
    METHOD          mp2
    BASIS           TODO_BASIS
    MEM_TOTAL       TODO_MEM
    UNRESTRICTED    hf_opt
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   scf_opt
    MAX_SCF_CYCLES  scf_cycle_opt
    N_FROZEN_CORE   fc
$end
```

The MP2 template is simpler — no CC convergence settings since MP2 is non-iterative
(it's a single perturbative correction on top of HF). The `N_FROZEN_CORE` variable
controls frozen core:

| Value | Meaning |
|-------|---------|
| `FC` (default) | Freeze the standard core orbitals (1s for C/N/O, 1s2s2p for Si/S) |
| `0` | Correlate all electrons (more expensive, rarely needed) |
| Integer `n` | Freeze exactly `n` lowest-energy orbitals |

> **Use the same `N_FROZEN_CORE` setting for both CCSD and MP2** runs on the same
> molecule so the amplitudes are directly comparable.

### 4.3 SLURM Generator (`cc_slurm_generator.sh` / `mp2_slurm_generator.sh`)

Both generators are identical shell scripts that loop over all `.in` files in the
current directory and create a `.slm` SLURM submission script for each:

```bash
for todo in *.in; do
    slm_file="${todo}.slm"
    cat > "$slm_file" <<EOL
#!/bin/bash
#SBATCH --job-name=$todo
#SBATCH --account=YOUR_ACCOUNT      ← your allocation
#SBATCH --time=                      ← wall-time limit (e.g., 04:00:00)
#SBATCH --nodes=                     ← number of nodes (usually 1)
#SBATCH -C cpu
#SBATCH --ntasks                     ← MPI tasks (usually 1 for Q-Chem)
#SBATCH --cpus-per-task              ← threads (e.g., 18 or 32)
#SBATCH -e ${todo}.err
#SBATCH --qos=regular

module load qchem
export QCLOCALSCR=\$QCSCRATCH

qchem $todo ${todo%.in}.out
EOL
done
```

**You need to fill in:**
- `--account`: your SLURM allocation/project code
- `--time`: estimated wall time (start generous, e.g., `08:00:00`)
- `--nodes`: `1` (Q-Chem runs on a single node for these calculations)
- `--ntasks`: `1`
- `--cpus-per-task`: number of OpenMP threads (match to `MEM_TOTAL` allocation)

---

## 5. End-to-End Pipeline

```
Step 1: Prepare Input Files
  ├── For each .xyz file:
  │   ├── Read charge & multiplicity (from comment line or table above)
  │   ├── Extract coordinates (lines 3+ of the XYZ file)
  │   ├── Choose method (CCSD or MP2) and basis set
  │   └── Fill in the template → save as molecule_name.in
  │
Step 2: Generate SLURM Scripts
  ├── Place all .in files in a directory
  ├── Run cc_slurm_generator.sh (or mp2_slurm_generator.sh)
  ├── Edit the generated .slm files with your account/resource details
  │
Step 3: Submit Jobs
  ├── sbatch molecule_name.in.slm
  │
Step 4: Monitor & Collect Output
  ├── Check .err files for failures
  ├── Check .out files for convergence
  │
Step 5: Parse Output
  ├── Extract energies (HF, MP2, CCSD, CCSD(T))
  ├── Extract t₁/t₂ amplitudes (leading amplitudes printed by default;
  │   full dump requires additional keywords — see Section 8)
  └── Extract MO coefficients (from .fchk file if GUI=2 is set)
```

---

## 6. Step-by-Step: Generating Input Files from XYZ

### Manual approach (good for understanding)

Take molecule `01_3_3_Sigmatropic_Rearrangements_...xyz`:

```
14                              ← atom count (ignore)
[3,3] Sigmatropic ...           ← comment line (ignore for input)
C          1.215804  -0.579377   0.252971    ← coordinates start here
C          1.471876   0.661540  -0.294775
O          0.357748  -1.382391  -0.255275
...
```

The Q-Chem input file becomes:

```
$molecule
0 1
C          1.215804  -0.579377   0.252971
C          1.471876   0.661540  -0.294775
O          0.357748  -1.382391  -0.255275
C         -1.330655  -0.780522   0.183694
C         -1.306470   0.533972  -0.295753
C         -0.474912   1.444998   0.328289
H          1.582783  -0.790571   1.266589
H          2.211156   1.306180   0.168621
H          1.252451   0.834759  -1.339760
H         -1.897208  -1.537806  -0.342408
H         -1.234077  -0.950402   1.250682
H         -1.612477   0.715924  -1.320269
H         -0.338659   2.440896  -0.078703
H         -0.279813   1.356482   1.390890
$end

$rem
    METHOD          ccsd(t)
    BASIS           STO-3G
    UNRESTRICTED    TRUE
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   DIIS_GDM
    MAX_SCF_CYCLES  200
    THRESH          15
    SCF_CONVERGENCE 10
    CC_CONVERGENCE  7
    GUI             2
$end
```

### Automation script

Here is a Python script to generate Q-Chem inputs from all your XYZ files:

```python
#!/usr/bin/env python3
"""Generate Q-Chem input files from XYZ coordinate files."""

import os
import glob
import argparse

CHARGE_MULT = {
    # ja4c17174
    "01_3_3_Sigmatropic": (0, 1),
    "02_Cationic_Aza_Cope": (1, 1),
    "03_TS_1_ketenimine_exo": (0, 1),
    "04_TS_SI_1_ketenimine_endo": (0, 1),
    "05_TS_2_formaldimine_5_exo": (0, 1),
    "06_TS_SI_2_formaldimine_5_endo": (0, 1),
    # ja025633n - all neutral singlets
}

CCSD_TEMPLATE = """$molecule
{charge} {mult}
{coords}
$end

$rem
    METHOD          {method}
    BASIS           {basis}
    {mem_line}UNRESTRICTED    {unrestricted}
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   DIIS_GDM
    MAX_SCF_CYCLES  200
    THRESH          15
    SCF_CONVERGENCE 10
    CC_CONVERGENCE  7
    GUI             2
$end
"""

MP2_TEMPLATE = """$molecule
{charge} {mult}
{coords}
$end

$rem
    METHOD          mp2
    BASIS           {basis}
    {mem_line}UNRESTRICTED    {unrestricted}
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   DIIS_GDM
    MAX_SCF_CYCLES  200
    N_FROZEN_CORE   FC
    GUI             2
$end
"""


def read_xyz(filepath):
    with open(filepath) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    coords = "".join(lines[2 : 2 + natoms])
    return natoms, coords.rstrip()


def get_charge_mult(filename):
    for key, (charge, mult) in CHARGE_MULT.items():
        if filename.startswith(key):
            return charge, mult
    return 0, 1  # default: neutral singlet


def generate_inputs(xyz_dir, output_dir, basis="STO-3G", method="ccsd(t)",
                    mem_mb=None, unrestricted=True):
    os.makedirs(output_dir, exist_ok=True)
    xyz_files = sorted(glob.glob(os.path.join(xyz_dir, "*.xyz")))

    for xyz_path in xyz_files:
        basename = os.path.splitext(os.path.basename(xyz_path))[0]
        natoms, coords = read_xyz(xyz_path)
        charge, mult = get_charge_mult(basename)

        mem_line = f"MEM_TOTAL       {mem_mb}\n    " if mem_mb else ""
        uhf = "TRUE" if unrestricted else "FALSE"

        if natoms > 40 or method.lower() == "mp2":
            template = MP2_TEMPLATE
            out_name = f"{basename}_mp2_{basis.replace('(', '').replace(')', '').replace('*', 's')}.in"
        else:
            template = CCSD_TEMPLATE
            out_name = f"{basename}_{basis.replace('(', '').replace(')', '').replace('*', 's')}.in"

        content = template.format(
            charge=charge, mult=mult, coords=coords,
            method=method, basis=basis, mem_line=mem_line,
            unrestricted=uhf,
        )

        out_path = os.path.join(output_dir, out_name)
        with open(out_path, "w") as f:
            f.write(content)
        print(f"  {out_name}  ({natoms} atoms, charge={charge}, mult={mult})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("xyz_dir", help="Directory containing .xyz files")
    parser.add_argument("output_dir", help="Directory for .in files")
    parser.add_argument("--basis", default="STO-3G")
    parser.add_argument("--method", default="ccsd(t)")
    parser.add_argument("--mem", type=int, default=None, help="MEM_TOTAL in MB")
    parser.add_argument("--restricted", action="store_true")
    args = parser.parse_args()

    generate_inputs(args.xyz_dir, args.output_dir, args.basis, args.method,
                    args.mem, not args.restricted)
```

**Usage:**

```bash
# Test run with minimal basis (STO-3G) — fast, validates the pipeline
python generate_qchem_inputs.py \
    coordinates_xyz_output/ja025633n_s1/ \
    inputs/ja025633n_sto3g/ \
    --basis STO-3G

# Production run with larger basis
python generate_qchem_inputs.py \
    coordinates_xyz_output/ja025633n_s1/ \
    inputs/ja025633n_631gd/ \
    --basis "6-31G*" --mem 24000

# MP2 for the large molecules
python generate_qchem_inputs.py \
    coordinates_xyz_output/ja4c17174_si_001/ \
    inputs/ja4c17174_mp2/ \
    --basis STO-3G --method mp2
```

---

## 7. Understanding the Output

The sample output (`sample_ccsd_output`) walks through a CCSD(T) calculation on
F₂ (fluorine dimer). Here's how to read each section:

### 7.1 SCF (Hartree-Fock) Convergence

```
 Cycle       Energy         DIIS error
    1    -198.5137076769      5.23e-02
    2    -198.6745719516      1.13e-02
    ...
   10    -198.6986752865      1.54e-09  Convergence criterion met
```

- The energy should decrease monotonically (roughly).
- DIIS error should drop by ~1 order of magnitude per cycle.
- **If it doesn't converge**: try `SCF_ALGORITHM = GDM`, increase
  `MAX_SCF_CYCLES`, or lower `SCF_CONVERGENCE`.

### 7.2 Orbital Information

```
Alpha MOs, Restricted
-- Occupied --
-26.447 -26.447  -1.775  -1.508  ...
-- Virtual --
  0.059   0.208   0.212   0.234  ...
```

- **Occupied orbitals**: filled with electrons. Negative energies (in Hartree).
- **Virtual orbitals**: empty, available for excitations. Positive energies.
- The gap between the highest occupied and lowest virtual orbital (HOMO-LUMO gap)
  indicates how "easy" the system is. Small gaps → harder SCF convergence, more
  multi-reference character.

### 7.3 Frozen Core

```
Alpha orbitals:
  - Frozen occupied    1    0    0    0    0    1    0    0     2
  - Active occupied    2    0    1    1    0    1    1    1     7
```

The 2 frozen orbitals are the 1s core electrons of each fluorine. Only the 7
active occupied orbitals participate in the correlation calculation.

### 7.4 CCSD Convergence

```
           Energy (a.u.)   Ediff      Tdiff       Comment
          -199.12680924                            ← initial MP2 guess
     1    -199.11826013   8.55e-03   5.34e-01
     ...
     9    -199.13419803   2.79e-07   7.12e-05
          -199.13419803                           CCSD T converged.
```

- **Ediff**: change in energy between iterations.
- **Tdiff**: change in t amplitudes between iterations.
- Convergence requires both to be below their thresholds.
- The initial guess comes from MP2 amplitudes (so MP2 is always computed as a
  byproduct of CCSD).

### 7.5 The Energy Summary (Most Important Part)

```
SCF energy                 = -198.69867529     ← Hartree-Fock energy
MP2 energy                 = -199.12680924     ← MP2 total energy
CCSD correlation energy    =   -0.43552274     ← correction on top of HF
CCSD total energy          = -199.13419803     ← HF + CCSD correlation
CCSD(T) correlation energy =   -0.01286438     ← perturbative triples
CCSD(T) total energy       = -199.14706241     ← HF + CCSD + (T)
```

Relationship: `E_CCSD(T) = E_HF + E_corr_CCSD + E_corr_(T)`

### 7.6 Leading Amplitudes

```
CCSD  T1^2 = 0.0055  T2^2 = 0.1179  Leading amplitudes:

Amplitude    Orbitals with energies
 0.0227       1 (B2u) A  ->  2 (B2u) A        ← single excitation (t₁)
             -0.8138           0.2340

-0.1294       3 (Ag) A   3 (Ag) B  ->  3 (B1u) A  3 (B1u) B   ← double excitation (t₂)
             -0.7569  -0.7569        0.0591   0.0591
```

- **T1^2**: sum of squares of all t₁ amplitudes. If T1^2 > 0.02, the system may
  have significant multi-reference character (HF is a poor starting point).
- **T2^2**: sum of squares of all t₂ amplitudes.
- Only the **leading** (largest) amplitudes are printed by default. To get the
  **full set** of amplitudes, additional keywords are needed (see next section).

---

## 8. Extracting Amplitudes and MO Coefficients

### 8.1 MO Coefficients via Formatted Checkpoint File

Add `GUI = 2` to the `$rem` section (already included in the templates above).
This generates a `.fchk` (formatted checkpoint) file containing:

- MO coefficients (the C matrix: AO → MO transformation)
- Orbital energies
- Overlap matrix
- Density matrices

The `.fchk` file is a plain-text file that can be parsed with standard tools
(PySCF, cclib, or custom scripts).

### 8.2 Full Amplitude Dump

The default output only shows **leading amplitudes**. To get the complete set of
t₁ and t₂ amplitudes, your collaborator mentioned:

> "Once you have qchem I will send where/how to modify to dump
> amplitudes/orbital coefficients."

This likely involves one or more of:
- `CC_PRINT` or related $rem variables to control verbosity of CC output
- Writing amplitudes to scratch files that can be read post-calculation
- Using the Q-Chem archive/qarchive format

**Action item**: Ask your collaborator for the specific keywords now that Q-Chem
is installed. Until then, the leading amplitudes from the default output and the
MO coefficients from the `.fchk` file give you a starting point.

### 8.3 Parsing the Output Programmatically

For energy extraction, a simple grep suffices:

```bash
# From a Q-Chem .out file:
grep "SCF energy"           molecule.out
grep "MP2 energy"           molecule.out
grep "CCSD total energy"    molecule.out
grep "CCSD(T) total energy" molecule.out
grep "T1\^2"                molecule.out
```

For the `.fchk` file, Python libraries like `cclib` or `pyscf` can parse it:

```python
import cclib
data = cclib.io.ccread("molecule.fchk")
mo_coeffs = data.mocoeffs      # shape: (nmo, nao)
mo_energies = data.moenergies  # in eV
```

---

## 9. Practical Tips and Troubleshooting

### SCF Won't Converge

1. **Try `SCF_ALGORITHM = GDM`** (more robust than DIIS, but slower).
2. **Increase `MAX_SCF_CYCLES`** to 500.
3. **Lower `SCF_CONVERGENCE`** from 10 to 8 (less tight).
4. **Add level shifting**: `LSHIFT = 200` (in units of 0.001 Hartree) stabilizes
   convergence by widening the HOMO-LUMO gap artificially.
5. **Check `INTERNAL_STABILITY = TRUE` output**: if the SCF solution is unstable,
   re-run with `UNRESTRICTED = TRUE` or try a broken-symmetry initial guess.

### Out of Memory

- For STO-3G: you shouldn't need `MEM_TOTAL` at all.
- For larger bases: set `MEM_TOTAL` to ~80% of your node's physical RAM in MB.
- If you hit memory limits: also add `MEM_STATIC = 5000` (reserves 5 GB for
  static allocations).
- The `(T)` triples correction is the most memory-hungry step. The output will
  warn you: "YOUR (T) CALCULATION MAY RUN OUT OF MEMORY".

### Which Molecules to Start With

Start with the **smallest molecules** and the **smallest basis set**:
1. `01_3_3_Sigmatropic` (14 atoms) with STO-3G — should finish in seconds.
2. `01_Unravelling` (16 atoms, ja025633n) with STO-3G.
3. Once confirmed working, scale up basis: STO-3G → 6-31G* → target.
4. Then tackle the 36-atom molecules.
5. Leave the 82-93 atom molecules for last (MP2 only).

### Basis Set Progression

| Stage | Basis | Purpose | Expected time (14 atoms) |
|-------|-------|---------|--------------------------|
| 0 | STO-3G | Pipeline validation | seconds |
| 1 | 6-31G(d) | Qualitative results | minutes |
| 2 | cc-pVDZ | Decent quantitative | minutes-hours |
| 3 | cc-pVTZ / def2-TZVPD | Production quality | hours-days |

---

## 10. Recommended Execution Plan

### Phase 1: Pipeline Validation (Day 1)

```bash
# 1. Generate STO-3G inputs for the small molecules
python generate_qchem_inputs.py \
    coordinates_xyz_output/ja4c17174_si_001/ \
    inputs/phase1_test/ \
    --basis STO-3G --method "ccsd(t)"

# 2. Run the 14-atom Claisen TS locally (no SLURM needed for STO-3G)
cd inputs/phase1_test/
qchem 01_3_3_Sigmatropic_STO-3G.in 01_3_3_Sigmatropic_STO-3G.out

# 3. Verify the output has all expected sections:
#    - SCF convergence
#    - CCSD convergence
#    - Energy summary
#    - Leading amplitudes
#    - .fchk file generated
```

### Phase 2: Small Molecules, Production Basis

Run CCSD(T) on the 14-16 atom molecules with 6-31G(d) and then the target basis.

### Phase 3: Medium Molecules (36 atoms)

Run CCSD(T)/STO-3G first to validate, then decide between CCSD(T) and MP2 for
larger bases based on resource availability.

### Phase 4: Large Molecules (82-93 atoms)

MP2 only. Even MP2 with a large basis on 93 atoms will be expensive. Discuss with
the collaborator about:
- Whether a smaller active space is acceptable
- Whether RI-MP2 (resolution-of-identity MP2) should be used for speedup
- Memory and time requirements

---

## 11. Questions to Clarify with Collaborator

Before the meeting, consider asking:

1. **Amplitude dump keywords**: You now have Q-Chem installed — what are the
   specific `$rem` variables or scratch file paths to dump the full t₁/t₂
   amplitude arrays (not just leading amplitudes)?

2. **Charge/multiplicity for molecules 03-06** (ja4c17174, the large Si-containing
   ones): Are these all neutral singlets, or do any carry a charge?

3. **CCSD vs CCSD(T)**: The template says `METHOD ccsd` but your checklist
   requests CCSD(T). Should we use `METHOD ccsd(t)` for all runs, or is CCSD
   sufficient for the amplitude initialization (since (T) doesn't produce
   additional amplitudes — it's a non-iterative energy correction)?

4. **MP2 amplitudes for LUCJ**: MP2 doesn't produce "t amplitudes" in the same
   iterative sense as CCSD. The MP2 "amplitudes" are the first-order doubles
   amplitudes t₂⁽¹⁾. Is there a specific Q-Chem keyword to dump these? Or do
   you extract them from the MP2 density matrix?

5. **RI-MP2 for large molecules**: Should we use resolution-of-identity MP2
   (`AUX_BASIS` keyword) for the 82-93 atom molecules to speed things up?

6. **Frozen core consistency**: The CCSD template doesn't explicitly set
   `N_FROZEN_CORE` (so Q-Chem uses the default `FC`). The MP2 template has `fc`
   as a placeholder. Should these be set to the same explicit value for
   consistency between CCSD and MP2 runs on the same molecule?

7. **Target basis set**: What is the final target basis? The checklist mentions
   both 6-31G(d) and def2-TZVPD. Is the plan to run both, or is one the primary
   target?

---

## Appendix A: Glossary

| Term | Definition |
|------|-----------|
| **SCF** | Self-Consistent Field — iterative procedure to find HF orbitals |
| **HF** | Hartree-Fock — mean-field approximation, no electron correlation |
| **MP2** | 2nd-order Møller-Plesset perturbation theory — cheapest correlation method |
| **CCSD** | Coupled Cluster with Singles and Doubles — systematic correlation method |
| **CCSD(T)** | CCSD + perturbative Triples — "gold standard" accuracy |
| **Basis set** | Set of functions used to expand molecular orbitals |
| **MO coefficients** | Matrix transforming atomic orbital basis into molecular orbitals |
| **t₁ amplitudes** | Single-excitation coefficients in coupled cluster |
| **t₂ amplitudes** | Double-excitation coefficients in coupled cluster |
| **Frozen core** | Excluding chemically inert inner-shell electrons from correlation |
| **DIIS** | Direct Inversion of Iterative Subspace — SCF acceleration algorithm |
| **GDM** | Geometric Direct Minimization — robust SCF algorithm for difficult cases |
| **LUCJ** | Locally-Unitary Coupled-Cluster Jastrow — quantum circuit ansatz |
| **Transition state (TS)** | Saddle point on the potential energy surface between reactants and products |
| **Single-point energy** | Energy calculation at a fixed geometry (no optimization) |
| **UHF** | Unrestricted Hartree-Fock — allows different spatial orbitals for alpha/beta spin |
| **`.fchk`** | Formatted checkpoint file — contains MO coefficients, energies, etc. |
| **Hartree (a.u.)** | Atomic unit of energy; 1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV |

## Appendix B: Quick Reference Card

```bash
# === VALIDATE PIPELINE ===
# Run smallest molecule with STO-3G (should take seconds)
qchem claisen_ts_sto3g.in claisen_ts_sto3g.out

# === CHECK OUTPUT ===
grep "Convergence criterion met" *.out    # SCF converged?
grep "CCSD T converged"         *.out    # CC converged?
grep "CCSD(T) total energy"     *.out    # Final energy
grep "T1\^2"                    *.out    # Multi-ref diagnostic

# === GENERATE SLURM JOBS ===
cd inputs/
bash ../templates/cc_slurm_generator.sh   # creates .slm for each .in
# Edit .slm files to fill in account, time, nodes, cpus
for f in *.slm; do sbatch "$f"; done      # submit all

# === MONITOR ===
squeue -u $USER                            # check job status
tail -20 molecule.out                      # peek at progress
```
