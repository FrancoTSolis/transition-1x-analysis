"""Command-line interface for the Q-Chem pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import ProjectConfig, SlurmParams, find_config


def _load_config(args) -> ProjectConfig:
    cfg_path = getattr(args, "config", None)
    if cfg_path:
        return ProjectConfig.from_yaml(cfg_path)
    found = find_config()
    if found:
        print(f"  Using config: {found}")
        return ProjectConfig.from_yaml(found)
    return ProjectConfig()


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_generate(args):
    """Generate Q-Chem input files from XYZ coordinates."""
    from .generate import generate_inputs

    config = _load_config(args)
    generate_inputs(
        args.xyz_dir,
        args.output_dir,
        config,
        method_override=args.method,
        basis_override=args.basis,
        mem_override=args.mem,
    )


def cmd_run(args):
    """Run Q-Chem jobs locally or via SLURM."""
    config = _load_config(args)

    if args.mode == "local":
        from .run import run_local_batch
        run_local_batch(
            args.input_dir,
            output_dir=args.output_dir,
            nprocs=args.nprocs,
            qchem_cmd=args.qchem_cmd,
            dry_run=args.dry_run,
        )
    elif args.mode == "slurm":
        slurm = config.slurm
        if args.account:
            slurm.account = args.account
        if args.time:
            slurm.time = args.time
        if args.nprocs:
            slurm.cpus_per_task = args.nprocs
        if not slurm.account:
            print("  Error: --account is required for SLURM mode (or set in pipeline_config.yaml)")
            sys.exit(1)

        from .run import submit_slurm_batch
        submit_slurm_batch(
            args.input_dir,
            slurm,
            output_dir=args.output_dir,
            nprocs=args.nprocs,
            dry_run=args.dry_run,
        )
    elif args.mode == "slurm-generate":
        slurm = config.slurm
        if args.account:
            slurm.account = args.account
        if args.time:
            slurm.time = args.time
        if args.nprocs:
            slurm.cpus_per_task = args.nprocs

        from .run import generate_slurm_scripts
        generate_slurm_scripts(
            args.input_dir,
            slurm,
            output_dir=args.output_dir,
            nprocs=args.nprocs,
        )


def cmd_parse(args):
    """Parse Q-Chem output files and report results."""
    from .parse import parse_batch, print_summary, export_json, export_csv

    results = parse_batch(args.output_dir)
    if not results:
        return

    print_summary(results)

    if args.json:
        export_json(results, args.json)
    if args.csv:
        export_csv(results, args.csv)


def cmd_pipeline(args):
    """Run the full pipeline: generate → run → parse."""
    from .generate import generate_inputs
    from .parse import parse_batch, print_summary, export_json, export_csv

    config = _load_config(args)

    # Step 1: Generate
    print("=" * 60)
    print("STEP 1: Generating Q-Chem input files")
    print("=" * 60)
    input_dir = Path(args.work_dir) / "inputs"
    generated = generate_inputs(
        args.xyz_dir,
        input_dir,
        config,
        method_override=args.method,
        basis_override=args.basis,
        mem_override=args.mem,
    )
    if not generated:
        print("  No inputs generated. Aborting.")
        sys.exit(1)

    # Step 2: Run
    print()
    print("=" * 60)
    print(f"STEP 2: Running jobs ({args.mode})")
    print("=" * 60)
    output_dir = Path(args.work_dir) / "outputs"

    if args.mode == "local":
        from .run import run_local_batch
        run_local_batch(
            input_dir,
            output_dir=output_dir,
            nprocs=args.nprocs,
            qchem_cmd=args.qchem_cmd,
            dry_run=args.dry_run,
        )
    elif args.mode == "slurm":
        slurm = config.slurm
        if args.account:
            slurm.account = args.account
        if args.time:
            slurm.time = args.time
        if not slurm.account:
            print("  Error: --account is required for SLURM mode")
            sys.exit(1)
        from .run import submit_slurm_batch
        submit_slurm_batch(input_dir, slurm, output_dir=output_dir, dry_run=args.dry_run)
        print("\n  Jobs submitted. Run `parse` after they complete.")
        return

    # Step 3: Parse (only for local mode)
    if not args.dry_run:
        print()
        print("=" * 60)
        print("STEP 3: Parsing results")
        print("=" * 60)
        results = parse_batch(output_dir)
        if results:
            print_summary(results)
            json_path = Path(args.work_dir) / "results.json"
            csv_path = Path(args.work_dir) / "results.csv"
            export_json(results, json_path)
            export_csv(results, csv_path)


def cmd_status(args):
    """Check status of SLURM jobs or scan output directory for results."""
    import subprocess
    import shutil

    out_dir = Path(args.output_dir)

    # Check for .out files and their completeness
    outs = sorted(out_dir.glob("*.out"))
    ins = sorted(out_dir.glob("*.in")) if not outs else []

    if outs:
        from .parse import parse_batch, print_summary
        results = parse_batch(out_dir)
        print_summary(results)
    else:
        # Check SLURM queue
        if shutil.which("squeue"):
            print("  No output files yet. Checking SLURM queue...\n")
            result = subprocess.run(
                ["squeue", "-u", str(subprocess.check_output(["whoami"]).decode().strip()),
                 "--format=%.18i %.30j %.8T %.10M %.6D %R"],
                capture_output=True, text=True,
            )
            if result.stdout.strip():
                print(result.stdout)
            else:
                print("  No SLURM jobs found in queue.")
        else:
            # Look for input files
            parent = out_dir.parent
            ins = sorted(parent.glob("inputs/*.in")) if (parent / "inputs").exists() else []
            if ins:
                print(f"  Found {len(ins)} input files, but no outputs yet.")
                print(f"  Input dir: {parent / 'inputs'}")
            else:
                print(f"  No .out files found in {out_dir}")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="qchem_pipeline",
        description="Q-Chem CCSD/MP2 pipeline for transition-state benchmarks",
    )
    parser.add_argument("--config", "-c", help="Path to pipeline_config.yaml")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- generate ---
    p_gen = sub.add_parser("generate", help="Generate Q-Chem input files from XYZ")
    p_gen.add_argument("xyz_dir", help="Directory containing .xyz files")
    p_gen.add_argument("output_dir", help="Directory for generated .in files")
    p_gen.add_argument("--method", "-m", help="Override method (e.g. ccsd(t), mp2)")
    p_gen.add_argument("--basis", "-b", help="Override basis set (e.g. STO-3G, 6-31G*)")
    p_gen.add_argument("--mem", type=int, help="MEM_TOTAL in MB")
    p_gen.set_defaults(func=cmd_generate)

    # --- run ---
    p_run = sub.add_parser("run", help="Execute Q-Chem jobs")
    p_run.add_argument("input_dir", help="Directory containing .in files")
    p_run.add_argument("--mode", choices=["local", "slurm", "slurm-generate"],
                       default="local",
                       help="Execution mode (default: local)")
    p_run.add_argument("--output-dir", "-o", help="Directory for output files (default: same as input)")
    p_run.add_argument("--nprocs", "-np", type=int, help="Number of CPU threads")
    p_run.add_argument("--qchem-cmd", default="qchem", help="Q-Chem executable (default: qchem)")
    p_run.add_argument("--account", help="SLURM account/allocation")
    p_run.add_argument("--time", help="SLURM wall time (e.g. 08:00:00)")
    p_run.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    p_run.set_defaults(func=cmd_run)

    # --- parse ---
    p_parse = sub.add_parser("parse", help="Parse Q-Chem output files")
    p_parse.add_argument("output_dir", help="Directory containing .out files")
    p_parse.add_argument("--json", help="Export results to JSON file")
    p_parse.add_argument("--csv", help="Export results to CSV file")
    p_parse.set_defaults(func=cmd_parse)

    # --- pipeline ---
    p_pipe = sub.add_parser("pipeline", help="Full pipeline: generate → run → parse")
    p_pipe.add_argument("xyz_dir", help="Directory containing .xyz files")
    p_pipe.add_argument("work_dir", help="Working directory for inputs/outputs/results")
    p_pipe.add_argument("--mode", choices=["local", "slurm"], default="local",
                        help="Execution mode (default: local)")
    p_pipe.add_argument("--method", "-m", help="Override method")
    p_pipe.add_argument("--basis", "-b", help="Override basis set")
    p_pipe.add_argument("--mem", type=int, help="MEM_TOTAL in MB")
    p_pipe.add_argument("--nprocs", "-np", type=int, help="Number of CPU threads")
    p_pipe.add_argument("--qchem-cmd", default="qchem", help="Q-Chem executable")
    p_pipe.add_argument("--account", help="SLURM account")
    p_pipe.add_argument("--time", help="SLURM wall time")
    p_pipe.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    p_pipe.set_defaults(func=cmd_pipeline)

    # --- status ---
    p_stat = sub.add_parser("status", help="Check job status and results")
    p_stat.add_argument("output_dir", help="Directory to check for .out files")
    p_stat.set_defaults(func=cmd_status)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
