#!/usr/bin/env python3
"""Move stale (spin-orbital) LUCJ outputs to a structured archive folder.

For each job directory, moves:
  lucj_*.npy, lucj_metadata.json, kappa_*.npy, .lucj_done
into stale_lucj/{job_name}/ while preserving the original CCSD outputs.

Usage:
    python archive_stale_lucj.py [--jobs-dir jobs] [--dry-run]
"""

import argparse
import os
import shutil
from pathlib import Path


LUCJ_FILES = [
    "lucj_diag_coulomb_mats.npy",
    "lucj_orbital_rotations.npy",
    "lucj_final_orbital_rotation.npy",
    "lucj_parameters.npy",
    "lucj_metadata.json",
    "kappa_real.npy",
    "kappa_imag.npy",
    ".lucj_done",
    ".lucj_failed",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs-dir", default="jobs")
    parser.add_argument("--archive-dir", default="stale_lucj_spin_orbital")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    jobs_dir = Path(args.jobs_dir)
    archive_dir = Path(args.archive_dir)

    moved = 0
    skipped = 0

    for d in sorted(jobs_dir.iterdir()):
        if not d.is_dir():
            continue

        has_lucj = any((d / f).exists() for f in LUCJ_FILES)
        if not has_lucj:
            skipped += 1
            continue

        dest = archive_dir / d.name
        if not args.dry_run:
            dest.mkdir(parents=True, exist_ok=True)

        for f in LUCJ_FILES:
            src = d / f
            if src.exists():
                if args.dry_run:
                    print(f"  {src} → {dest / f}")
                else:
                    shutil.move(str(src), str(dest / f))

        moved += 1

    print(f"{'Would move' if args.dry_run else 'Moved'} LUCJ files from {moved} jobs")
    print(f"Skipped {skipped} jobs (no LUCJ files)")
    if not args.dry_run:
        print(f"Archive location: {archive_dir}/")


if __name__ == "__main__":
    main()
