#!/usr/bin/env python3
"""Precompute kappa = logm(U) for all job directories.

Saves kappa_real.npy and kappa_imag.npy in each job directory.
Runs in parallel to speed up the 30k jobs.

Usage:
    python -m pretrain.precompute_kappa --jobs-dir jobs --workers 12
"""

import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.linalg import logm


def process_job(job_dir: str) -> str:
    job_dir = Path(job_dir)
    name = job_dir.name

    kr_path = job_dir / "kappa_real.npy"
    ki_path = job_dir / "kappa_imag.npy"

    if kr_path.exists() and ki_path.exists():
        return f"SKIP {name}"

    U_path = job_dir / "lucj_orbital_rotations.npy"
    if not U_path.exists():
        return f"MISS {name}"

    try:
        U = np.load(U_path)
        kappa_real = np.zeros_like(U.real)
        kappa_imag = np.zeros_like(U.imag)
        for rep in range(U.shape[0]):
            kappa = logm(U[rep])
            kappa_real[rep] = kappa.real
            kappa_imag[rep] = kappa.imag

        np.save(kr_path, kappa_real.astype(np.float64))
        np.save(ki_path, kappa_imag.astype(np.float64))
        return f"DONE {name}"
    except Exception as e:
        return f"FAIL {name}: {e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs-dir", default="jobs")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()

    jobs_dir = Path(args.jobs_dir)
    job_dirs = sorted(
        str(d) for d in jobs_dir.iterdir()
        if d.is_dir() and (d / ".lucj_done").exists()
    )
    print(f"Found {len(job_dirs)} completed LUCJ jobs")

    done = skip = fail = 0
    with Pool(args.workers) as pool:
        for i, result in enumerate(pool.imap_unordered(process_job, job_dirs)):
            if result.startswith("DONE"):
                done += 1
            elif result.startswith("SKIP"):
                skip += 1
            else:
                fail += 1
                print(result)

            if (i + 1) % 1000 == 0:
                print(f"  {i+1}/{len(job_dirs)}: done={done}, skip={skip}, fail={fail}")

    print(f"\nFinished: done={done}, skip={skip}, fail={fail}")


if __name__ == "__main__":
    main()
