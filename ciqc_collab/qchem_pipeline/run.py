"""Execute Q-Chem calculations locally or via SLURM."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .config import SlurmParams


# ---------------------------------------------------------------------------
# Local execution
# ---------------------------------------------------------------------------

def run_local(
    input_file: str | Path,
    output_dir: str | Path | None = None,
    *,
    nprocs: int | None = None,
    qchem_cmd: str = "qchem",
    dry_run: bool = False,
) -> Path:
    """Run a single Q-Chem job directly via subprocess.

    Returns the path to the output file.
    """
    inp = Path(input_file).resolve()
    if not inp.exists():
        raise FileNotFoundError(inp)

    out_dir = Path(output_dir).resolve() if output_dir else inp.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / inp.with_suffix(".out").name

    cmd = [qchem_cmd]
    if nprocs:
        cmd += ["-nt", str(nprocs)]
    cmd += [str(inp), str(out_file)]

    if dry_run:
        print(f"  [dry-run] {' '.join(cmd)}")
        return out_file

    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(out_dir))
    if result.returncode != 0:
        err_msg = result.stderr[:2000] if result.stderr else "(no stderr)"
        print(f"  FAILED (exit {result.returncode}): {err_msg}")
    else:
        print(f"  Done → {out_file}")

    # Q-Chem may drop .fchk in cwd or next to input; move stray copies to out_dir
    for candidate_dir in [inp.parent, Path.cwd()]:
        fchk = candidate_dir / inp.with_suffix(".fchk").name
        if fchk.exists() and fchk.parent != out_dir:
            dest = out_dir / fchk.name
            shutil.move(str(fchk), str(dest))
            print(f"  Moved {fchk.name} → {out_dir}")

    return out_file


def run_local_batch(
    input_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    nprocs: int | None = None,
    qchem_cmd: str = "qchem",
    dry_run: bool = False,
) -> list[Path]:
    """Run all .in files in *input_dir* sequentially."""
    input_dir = Path(input_dir)
    inputs = sorted(input_dir.glob("*.in"))
    if not inputs:
        print(f"  No .in files found in {input_dir}")
        return []
    print(f"  Found {len(inputs)} input files in {input_dir}\n")
    outputs = []
    for inp in inputs:
        out = run_local(inp, output_dir, nprocs=nprocs, qchem_cmd=qchem_cmd, dry_run=dry_run)
        outputs.append(out)
    return outputs


# ---------------------------------------------------------------------------
# SLURM submission
# ---------------------------------------------------------------------------

def _build_slurm_script(
    input_file: Path,
    output_file: Path,
    params: SlurmParams,
    *,
    nprocs: int | None = None,
    run_dir: Path | None = None,
) -> str:
    """Build a SLURM batch script for a single Q-Chem job.

    *run_dir* is the directory the script will ``cd`` into before running
    Q-Chem, so that ``.fchk`` files land alongside the ``.out`` files.
    """
    job_name = input_file.stem
    cpus = nprocs or params.cpus_per_task

    header_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
    ]
    if params.account:
        header_lines.append(f"#SBATCH --account={params.account}")
    header_lines += [
        f"#SBATCH --time={params.time}",
        f"#SBATCH --nodes={params.nodes}",
        f"#SBATCH --ntasks={params.ntasks}",
        f"#SBATCH --cpus-per-task={cpus}",
    ]
    if params.constraint:
        header_lines.append(f"#SBATCH -C {params.constraint}")
    header_lines += [
        f"#SBATCH --qos={params.qos}",
        f"#SBATCH -o {job_name}.slurm.log",
        f"#SBATCH -e {job_name}.slurm.err",
    ]

    body_lines = [""]
    for mod in params.modules:
        body_lines.append(f"module load {mod}")
    body_lines.append("export QCLOCALSCR=${QCSCRATCH:-/tmp}")
    for k, v in params.extra_exports.items():
        body_lines.append(f"export {k}={v}")

    if run_dir:
        body_lines += ["", f"cd {run_dir}"]

    body_lines += [
        "",
        f"qchem -nt {cpus} {input_file.resolve()} {output_file.name}",
        "",
        f'echo "Finished {job_name} at $(date)"',
    ]

    return "\n".join(header_lines + body_lines) + "\n"


def generate_slurm_scripts(
    input_dir: str | Path,
    params: SlurmParams,
    *,
    output_dir: str | Path | None = None,
    nprocs: int | None = None,
) -> list[Path]:
    """Create a .slm SLURM script for each .in file.

    Returns paths to generated scripts.
    """
    input_dir = Path(input_dir)
    out_dir = Path(output_dir) if output_dir else input_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    inputs = sorted(input_dir.glob("*.in"))
    if not inputs:
        print(f"  No .in files found in {input_dir}")
        return []

    scripts: list[Path] = []
    for inp in inputs:
        out_file = out_dir / (inp.stem + ".out")
        content = _build_slurm_script(
            inp, out_file, params, nprocs=nprocs, run_dir=out_dir.resolve(),
        )
        slm_path = out_dir / (inp.stem + ".slm")
        slm_path.write_text(content)
        scripts.append(slm_path)
        print(f"  {slm_path.name}")

    print(f"\n  Generated {len(scripts)} SLURM scripts in {out_dir}")
    return scripts


def submit_slurm_batch(
    input_dir: str | Path,
    params: SlurmParams,
    *,
    output_dir: str | Path | None = None,
    nprocs: int | None = None,
    dry_run: bool = False,
) -> list[str]:
    """Generate SLURM scripts and submit them via sbatch.

    Returns a list of SLURM job IDs.
    """
    if not shutil.which("sbatch"):
        raise RuntimeError("sbatch not found — are you on a SLURM cluster?")

    scripts = generate_slurm_scripts(input_dir, params, output_dir=output_dir, nprocs=nprocs)
    job_ids: list[str] = []

    for slm in scripts:
        if dry_run:
            print(f"  [dry-run] sbatch {slm}")
            job_ids.append("DRY_RUN")
            continue

        result = subprocess.run(
            ["sbatch", str(slm)],
            capture_output=True, text=True,
            cwd=str(slm.parent),
        )
        if result.returncode == 0:
            # sbatch output: "Submitted batch job 12345"
            parts = result.stdout.strip().split()
            jid = parts[-1] if parts else "?"
            job_ids.append(jid)
            print(f"  Submitted {slm.name} → job {jid}")
        else:
            print(f"  FAILED to submit {slm.name}: {result.stderr.strip()}")
            job_ids.append("FAILED")

    return job_ids
