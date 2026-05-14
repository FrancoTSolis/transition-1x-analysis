# CCSD Amplitude Extraction -- Transition1x Dataset

Extract full CCSD T1/T2 amplitude tensors for **reactant**, **transition state**,
and **product** geometries from the Transition1x dataset using a custom-patched
Q-Chem (with `PRINT_QIS` support).

## Overview

The Transition1x dataset contains 10,073 reactions, each with dedicated
single-frame geometries for:
- **Reactant** (R)
- **Transition State** (TS)
- **Product** (P)

This subfolder generates Q-Chem CCSD/STO-3G inputs for all three species of
every reaction and extracts the full amplitude tensors needed to initialize
LUCJ (Locally-Unitary Coupled-Cluster Jastrow) ansatz on a quantum circuit.

### Dataset Statistics

| Property | Value |
|----------|-------|
| Total reactions | 10,073 |
| Total Q-Chem jobs | 30,219 (3 per reaction) |
| Atom count range | 4--23 |
| Basis functions (STO-3G) | 18--51 |
| Occupied orbitals (avg) | ~18 |
| Virtual orbitals (avg) | ~18 |

### Disk Space Estimates

| Component | Size |
|-----------|------|
| Amplitude files (T1 + T2) | ~106 GB |
| Q-Chem outputs (.out + .fchk) | ~5 GB |
| XYZ + input files | ~1 GB |
| **Total permanent storage** | **~112 GB** |
| Scratch (per concurrent job) | ~100--400 MB |
| Peak scratch (48 concurrent) | ~5--20 GB |
| Available disk space | **350 TB** |

**Disk space is not a concern** on this filesystem (350 TB available, 4% used).

### Timing Estimates

Benchmarked on scai2 (48 cores, 2x Intel Xeon E5-2650 v4 @ 2.20 GHz, 220 GB RAM):

| Molecule size | Wall time (1 thread) | Scratch size |
|--------------|---------------------|-------------|
| 7 atoms (27 BFs) | ~27s | ~100 MB |
| 14 atoms (36 BFs) | ~2 min | ~200 MB |
| 18 atoms (46 BFs) | ~3.3 min | ~222 MB |
| 23 atoms (51 BFs) | ~10--30 min | ~400 MB |

**Key finding**: For STO-3G calculations, parallelism within a single Q-Chem
job provides **no speedup** -- the problem is too small for OpenMP threading
overhead to be worthwhile. The optimal strategy is:

> **Run many single-threaded jobs in parallel.**

With 48 cores running 48 single-threaded jobs concurrently, estimated
throughput is ~6,500 jobs/hour for the smallest molecules, scaling down
to ~200 jobs/hour for the largest. Weighted average across the dataset
is approximately **2,000--3,000 jobs/hour**, meaning the full 30,219 jobs
should complete in **~10--15 hours**.

## Prerequisites

1. **Custom Q-Chem** with `PRINT_QIS` amplitude dump patch:
   ```bash
   source /xuanwu-tank/east/fts/qchem_compile/qcenv_dev.sh
   ```

2. **Python 3** with `h5py` and `numpy` (for geometry extraction):
   ```bash
   pip install h5py numpy
   ```

3. **Transition1x HDF5 file** at `../data/transition1x.h5`

## Step-by-Step Instructions

### 1. Generate Q-Chem Input Files

```bash
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes

# Preview what will be generated (no files written)
python3 extract_geometries.py --dry-run

# Generate all 30,219 jobs
python3 extract_geometries.py

# Or generate a subset (e.g., molecules with <=14 atoms)
python3 extract_geometries.py --max-atoms 14

# Or a single formula / reaction
python3 extract_geometries.py --formula C2H2N2O
python3 extract_geometries.py --formula C2H2N2O --rxn-id rxn2091
```

Options: `--h5-path`, `--output-dir` (default: `jobs`), `--max-atoms`,
`--min-atoms`, `--max-reactions`, `--split {data,train,test,val}`, `--dry-run`.

### 2. Run a Single Job (Test)

```bash
bash run_single.sh jobs/C2H2N2O_rxn2091_TS 1
```

This runs Q-Chem CCSD with `PRINT_QIS` + `-save`, copies amplitude `.dat`
files from scratch to the job directory, cleans up scratch, and creates a
`.done` or `.failed` marker.

### 3. Run All Jobs in Parallel (Background)

**Recommended** (full speed, 48-core machine):

```bash
nohup bash run_batch.sh \
    --parallel-jobs 48 \
    --threads-per-job 1 \
    > batch_run.log 2>&1 &
disown
```

**Conservative** (leaves cores for other work):

```bash
nohup bash run_batch.sh \
    --parallel-jobs 24 \
    --threads-per-job 1 \
    > batch_run.log 2>&1 &
disown
```

**Test run** (first 100 jobs only):

```bash
bash run_batch.sh --max-jobs 100 --parallel-jobs 24 --threads-per-job 1
```

**Dry run** (preview without executing):

```bash
bash run_batch.sh --dry-run
```

#### Batch Runner Options

| Option | Default | Description |
|--------|---------|-------------|
| `--jobs-dir DIR` | `jobs` | Directory containing job subdirs |
| `--threads-per-job N` | 4 | OpenMP threads per Q-Chem process |
| `--parallel-jobs N` | 6 | Number of simultaneous jobs |
| `--log FILE` | `run_batch.log` | Log file path |
| `--max-jobs N` | 0 (all) | Limit jobs to process |
| `--dry-run` | -- | Print job list, don't execute |

### 4. Monitor Progress

```bash
bash check_progress.sh          # Quick summary
tail -f run_batch.log           # Watch log in real time
find jobs -name ".done" | wc -l # Count completed
```

### 5. Resume After Interruption

The batch runner is **idempotent** -- it skips jobs with a `.done` marker.
Just re-run `run_batch.sh`. Failed jobs (`.failed` marker) are auto-retried.

### 6. Handle Failed Jobs

```bash
find jobs -name ".failed" -exec dirname {} \;   # List failures
find jobs -name ".failed" -delete               # Clear markers
bash run_batch.sh --parallel-jobs 48 --threads-per-job 1  # Retry
```

## Output File Format

### Amplitude Files

Each `.dat` file is plain text:

```
T1 amplitudes:
2 18 18
   5.9865477017612861e-04    1.6743251334230430e-03  ...
```

- **Line 1**: Header (`T1 amplitudes:`, `T2 amplitudes:`, `MP2 T2 amplitudes:`)
- **Line 2**: Tensor rank + dimensions (e.g., `2 18 18` = rank-2, 18x18)
- **Remaining**: Elements in row-major order, 3 per line, full double precision

### Manifest

`jobs/manifest.json` contains metadata for every reaction including formula,
rxn_id, natoms, n_electrons, n_basis, n_occ, n_virt, and reference energies
(wB97x/6-31G(d) in eV) for all three species.

## Parallelism Benchmark Results

Tested on 7-atom molecule (C2H2N2O, 27 basis functions, STO-3G):

| Threads | Wall (s) | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1 | 26.5 | 1.00x | 100% |
| 2 | 71.8 | 0.37x | 18% |
| 4 | 58.9 | 0.45x | 11% |
| 8 | 65.9 | 0.40x | 5% |
| 12 | 66.9 | 0.40x | 3% |
| 16 | 45.9 | 0.58x | 4% |
| 24 | 46.5 | 0.57x | 2% |

Multi-threading actually **hurts** performance for these small STO-3G molecules.
Startup and synchronization overhead of OpenMP exceeds any computational
benefit. **Optimal strategy: 1 thread per job, maximize parallel jobs.**

| Configuration | Cores used | Throughput |
|--------------|-----------|-----------|
| 1t x 48 parallel | 48 | ~6,500 jobs/hr |
| 2t x 24 parallel | 48 | ~1,200 jobs/hr |
| 4t x 12 parallel | 48 | ~730 jobs/hr |

## Transition1x Data Structure

Each reaction in the HDF5 has dedicated single-geometry entries -- no need
to parse trajectories:

```
data/{formula}/{rxn_id}/
    atomic_numbers          # shape (n_atoms,)
    positions               # shape (n_frames, n_atoms, 3) -- full trajectory
    reactant/               # <-- single-frame dedicated geometry
        positions           # shape (1, n_atoms, 3)
        wB97x_6-31G(d).energy
    transition_state/       # <-- single-frame dedicated geometry
    product/                # <-- single-frame dedicated geometry
```

The `positions` at the reaction level is the full trajectory (hundreds of
frames). The `reactant`, `transition_state`, and `product` subgroups each
contain a **single optimized geometry** -- these are what we extract.

## NFS Performance: Why `run_single.sh` Uses `cd`

Q-Chem's launch scripts (`serial.csh`, `qchem`) create temporary files
(`.qcin`, `.out.files/`, `.in.fchk`) in the **current working directory**.
If Q-Chem is invoked from a directory containing tens of thousands of entries
(like `jobs/` with 30k subdirectories), every `open()` / `stat()` / `creat()`
triggers an NFS LOOKUP against that huge directory. With multiple processes
doing this concurrently, NFS serializes the lookups and throughput collapses.

In testing, running 12 parallel jobs from the top-level `ccsd_amplitudes/`
directory resulted in **~390s per job**. After changing `run_single.sh` to
`cd` into each job's own directory (which contains only 2--3 files) before
invoking Q-Chem, the same jobs completed in **~10--25s** -- a 15--40x
improvement. This also keeps all Q-Chem artifacts (`.fchk`, `.out.files/`,
etc.) contained inside each job directory rather than littering the parent.

## LUCJ Compressed Double Factorization (ffsim)

After the CCSD amplitude extraction, a second stage computes **LUCJ
(Locally-Unitary Coupled-Cluster Jastrow) operators** from the T1/T2
amplitudes using ffsim's compressed double factorization. This produces
the U (orbital rotation) and J (diagonal Coulomb) matrices needed to
initialize a quantum circuit ansatz.

### What It Computes

For each job directory that contains `ccsd_t1.dat` and `ccsd_t2.dat`,
`compute_lucj.py` calls:

```python
lucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
    t2=ccsd_t2,
    t1=ccsd_t1,
    n_reps=2,
    interaction_pairs=(pairs_aa, pairs_ab),
    optimize=True,
    options=dict(maxiter=50),
)
```

where `pairs_aa = [(p, p+1) for p in range(norb-1)]` and
`pairs_ab = [(p, p) for p in range(norb)]` model square-lattice hardware
connectivity. `n_reps=2` gives two circuit layers (as required by QSD).
`optimize=True` enables the compressed double factorization, which adjusts
the truncated U/J matrices to better recover the original T2 amplitudes.

### Output Files (Per Job Directory)

| File | Description |
|------|-------------|
| `lucj_diag_coulomb_mats.npy` | J matrices, shape `(n_reps, 2, norb, norb)` |
| `lucj_orbital_rotations.npy` | U matrices, shape `(n_reps, norb, norb)` |
| `lucj_final_orbital_rotation.npy` | Final orbital rotation, shape `(norb, norb)` |
| `lucj_parameters.npy` | Flat parameter vector for the LUCJ operator |
| `lucj_metadata.json` | Metadata: norb, n_reps, nocc, nvirt, timing, etc. |
| `.lucj_done` | Completion marker |

### Running

```bash
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes

# Single job (must run from outside the project root to avoid
# the local ffsim/ source tree shadowing the pip-installed package):
cd /tmp
python3 /path/to/compute_lucj.py jobs/C2H2N2O_rxn2091_TS

# Full batch (recommended: 6 parallel, 6 threads each = 36 of 48 cores):
nohup bash run_lucj_batch.sh --parallel-jobs 6 --threads-per-job 6 \
    > lucj_batch.log 2>&1 &
disown

# Monitor progress:
grep -c "DONE" run_lucj_batch.log
grep "DONE" run_lucj_batch.log | tail -5
```

### Parallelism

Unlike the Q-Chem CCSD stage (which is single-threaded and I/O-bound),
the ffsim LUCJ computation is **CPU-bound and internally multithreaded**
via jax/numpy. Each job uses ~6 CPU cores for ~30--55s per molecule
(depending on orbital count). The optimal configuration on this 48-core
machine is:

| Configuration | Cores used | Per-job wall time |
|--------------|-----------|-------------------|
| 6 parallel x 6 threads | 36 | ~28--55s |

The batch runner sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`MKL_NUM_THREADS`, and `JAX_NUM_THREADS` to control thread usage per job.

### Note on ffsim Import

The project root contains a local `ffsim/` source checkout (used for
reference and the tutorial notebook). This **shadows** the pip-installed
ffsim package. The batch runner avoids this by `cd /tmp` before invoking
Python. If running `compute_lucj.py` manually, do the same.

## Notes

- The custom Q-Chem dev build (SVN trunk r47985) does not require a license.
- `-save` flag is **mandatory** to preserve the scratch directory where
  amplitude files are written.
- Scratch is cleaned up after each job to avoid filling the scratch partition.
- The `.done` / `.failed` marker system makes runs fully resumable.
- Energies from Transition1x are in **eV** (wB97x/6-31G(d)); CCSD energies
  in our outputs are in **Hartree** (1 Ha = 27.211386 eV).
