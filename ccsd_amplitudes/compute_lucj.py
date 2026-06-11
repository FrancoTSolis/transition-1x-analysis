#!/usr/bin/env python3
"""Compute LUCJ U and J matrices from CCSD t-amplitudes using ffsim.

Reads ccsd_t1.dat and ccsd_t2.dat from a job directory, computes the
compressed double-factorized LUCJ operator via ffsim, and saves the
resulting diag_coulomb_mats (J), orbital_rotations (U),
final_orbital_rotation, and flat parameter vector.

All I/O is confined to the job directory to avoid race conditions
when running many jobs in parallel.

Usage:
    python compute_lucj.py <job_dir>
    python compute_lucj.py <job_dir> --n-reps 2 --maxiter 50
    python compute_lucj.py --help

IMPORTANT: This script must NOT be run from a directory that contains
a local 'ffsim' folder (e.g. the project root), as that would shadow
the installed ffsim package. The batch runner handles this automatically.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


def load_amplitude_file(path: Path) -> np.ndarray:
    """Load a Q-Chem amplitude .dat file into a numpy array."""
    lines = path.read_text().splitlines()
    header = lines[0].strip()
    shape_line = lines[1].strip().split()
    rank = int(shape_line[0])
    dims = tuple(int(d) for d in shape_line[1:])
    assert len(dims) == rank, f"Rank mismatch: {rank} vs {len(dims)} in {path}"

    values = []
    for line in lines[2:]:
        for token in line.strip().split():
            values.append(float(token))

    total = 1
    for d in dims:
        total *= d
    assert len(values) == total, (
        f"Element count mismatch in {path}: expected {total}, got {len(values)}"
    )

    return np.array(values).reshape(dims)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("job_dir", help="Path to job directory containing .dat files")
    parser.add_argument("--n-reps", type=int, default=2,
                        help="Number of LUCJ repetitions (circuit layers, default: 2)")
    parser.add_argument("--maxiter", type=int, default=50,
                        help="Max iterations for compressed DF optimization (default: 50)")
    parser.add_argument("--regularization", type=float, default=0.0,
                        help="Regularization parameter for compressed DF (default: 0.0)")
    args = parser.parse_args()

    job_dir = Path(args.job_dir).resolve()
    if not job_dir.is_dir():
        print(f"ERROR: {job_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    done_marker = job_dir / ".lucj_done"
    fail_marker = job_dir / ".lucj_failed"

    if done_marker.exists():
        print(f"SKIP: {job_dir.name} LUCJ already computed")
        return

    fail_marker.unlink(missing_ok=True)

    t1_path = job_dir / "ccsd_t1.dat"
    t2_path = job_dir / "ccsd_t2.dat"

    if not t1_path.exists() or not t2_path.exists():
        msg = f"SKIP: {job_dir.name} missing amplitude files"
        print(msg)
        fail_marker.write_text("MISSING_AMPLITUDES\n")
        return

    try:
        import ffsim
    except ImportError:
        print("ERROR: ffsim not installed or not importable", file=sys.stderr)
        sys.exit(1)

    start = time.time()
    print(f"START {job_dir.name} (n_reps={args.n_reps}, maxiter={args.maxiter})")

    t1_spinorb = load_amplitude_file(t1_path)
    t2_spinorb = load_amplitude_file(t2_path)

    nocc_spin, nvirt_spin = t1_spinorb.shape
    nocc = nocc_spin // 2
    nvirt = nvirt_spin // 2
    norb = nocc + nvirt

    t1 = t1_spinorb[:nocc, :nvirt]
    t2 = t2_spinorb[:nocc, :nocc, :nvirt, :nvirt]

    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    pairs_ab = [(p, p) for p in range(norb)]

    try:
        lucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
            t2=t2,
            t1=t1,
            n_reps=args.n_reps,
            interaction_pairs=(pairs_aa, pairs_ab),
            optimize=True,
            options=dict(maxiter=args.maxiter),
            regularization=args.regularization,
        )
    except Exception as e:
        elapsed = time.time() - start
        msg = f"FFSIM_ERROR: {e}"
        print(f"FAIL  {job_dir.name}  {msg}  wall={elapsed:.1f}s")
        fail_marker.write_text(f"{msg}\n")
        return

    np.save(job_dir / "lucj_diag_coulomb_mats.npy", lucj_op.diag_coulomb_mats)
    np.save(job_dir / "lucj_orbital_rotations.npy", lucj_op.orbital_rotations)
    if lucj_op.final_orbital_rotation is not None:
        np.save(job_dir / "lucj_final_orbital_rotation.npy",
                lucj_op.final_orbital_rotation)
    np.save(job_dir / "lucj_parameters.npy", lucj_op.to_parameters())

    metadata = {
        "norb": int(lucj_op.norb),
        "n_reps": int(lucj_op.n_reps),
        "nocc": int(nocc),
        "nvirt": int(nvirt),
        "nocc_spin": int(nocc_spin),
        "nvirt_spin": int(nvirt_spin),
        "basis": "spatial",
        "n_params": int(len(lucj_op.to_parameters())),
        "diag_coulomb_mats_shape": list(lucj_op.diag_coulomb_mats.shape),
        "orbital_rotations_shape": list(lucj_op.orbital_rotations.shape),
        "maxiter": args.maxiter,
        "regularization": args.regularization,
        "wall_time_s": round(time.time() - start, 2),
    }
    (job_dir / "lucj_metadata.json").write_text(json.dumps(metadata, indent=2))

    done_marker.touch()
    elapsed = time.time() - start
    print(f"DONE  {job_dir.name}  norb={norb} n_reps={args.n_reps} "
          f"n_params={metadata['n_params']}  wall={elapsed:.1f}s")


if __name__ == "__main__":
    main()
