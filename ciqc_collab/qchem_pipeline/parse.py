"""Parse Q-Chem output files for energies, convergence, and amplitudes."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field, asdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class SCFResult:
    converged: bool = False
    energy: float | None = None
    cycles: int = 0
    wall_time: str = ""


@dataclass
class CCResult:
    converged: bool = False
    ccsd_energy: float | None = None
    ccsd_corr_energy: float | None = None
    ccsd_t_energy: float | None = None
    ccsd_t_corr_energy: float | None = None
    t1_squared: float | None = None
    t2_squared: float | None = None
    iterations: int = 0


@dataclass
class MP2Result:
    energy: float | None = None


@dataclass
class Amplitude:
    """A single leading amplitude from the output."""
    value: float
    excitation_type: str  # "single" or "double"
    from_orbitals: str
    to_orbitals: str


@dataclass
class ParsedOutput:
    filename: str
    success: bool = False
    error: str | None = None
    scf: SCFResult = field(default_factory=SCFResult)
    mp2: MP2Result = field(default_factory=MP2Result)
    cc: CCResult = field(default_factory=CCResult)
    leading_amplitudes: list[Amplitude] = field(default_factory=list)
    fchk_exists: bool = False
    method: str = ""
    basis: str = ""
    natoms: int = 0
    nbasis: int = 0
    wall_time: str = ""


# ---------------------------------------------------------------------------
# Parsing engine
# ---------------------------------------------------------------------------

_RE_SCF_ENERGY = re.compile(r"SCF\s+energy\s*=\s*([-\d.]+)")
_RE_MP2_ENERGY = re.compile(r"MP2\s+energy\s*=\s*([-\d.]+)")
_RE_CCSD_TOTAL = re.compile(r"CCSD\s+total\s+energy\s*=\s*([-\d.]+)")
_RE_CCSD_CORR = re.compile(r"CCSD\s+correlation\s+energy\s*=\s*([-\d.]+)")
_RE_CCSD_T_TOTAL = re.compile(r"CCSD\(T\)\s+total\s+energy\s*=\s*([-\d.]+)")
_RE_CCSD_T_CORR = re.compile(r"CCSD\(T\)\s+correlation\s+energy\s*=\s*([-\d.]+)")
_RE_T_DIAG = re.compile(r"T1\^2\s*=\s*([\d.]+)\s+T2\^2\s*=\s*([\d.]+)")
_RE_SCF_CONVERGED = re.compile(r"Convergence criterion met")
_RE_CC_CONVERGED = re.compile(r"CCSD T converged")
_RE_SCF_CYCLE = re.compile(r"^\s+(\d+)\s+([-\d.]+)\s+([\d.eE+-]+)", re.MULTILINE)
_RE_BASIS_FUNCS = re.compile(r"There are\s+\d+\s+shells and\s+(\d+)\s+basis functions")
_RE_NATOMS_LINE = re.compile(r"Standard Nuclear Orientation.*\n.*I\s+Atom.*\n.*-+\n((?:.*\n)*?).*-+")
_RE_TOTAL_TIME = re.compile(r"Total job time:\s*([\d.]+s\(wall\))")
_RE_METHOD = re.compile(r"METHOD\s+(\S+)", re.IGNORECASE)
_RE_BASIS = re.compile(r"BASIS\s+(\S+)", re.IGNORECASE)
_RE_AMPLITUDE_SINGLE = re.compile(
    r"\s+([-\d.]+)\s+(\d+\s*\([^)]+\)\s*[AB])\s*->\s*(\d+\s*\([^)]+\)\s*[AB])"
)
_RE_AMPLITUDE_DOUBLE = re.compile(
    r"\s+([-\d.]+)\s+(\d+\s*\([^)]+\)\s*[AB]\s+\d+\s*\([^)]+\)\s*[AB])\s*->\s*(\d+\s*\([^)]+\)\s*[AB]\s+\d+\s*\([^)]+\)\s*[AB])"
)


def parse_output(path: str | Path) -> ParsedOutput:
    """Parse a Q-Chem .out file and extract key results."""
    path = Path(path)
    result = ParsedOutput(filename=path.name)

    if not path.exists():
        result.error = "File not found"
        return result

    text = path.read_text(errors="replace")
    if not text.strip():
        result.error = "Empty file"
        return result

    # Method & basis
    m = _RE_METHOD.search(text)
    if m:
        result.method = m.group(1)
    m = _RE_BASIS.search(text)
    if m:
        result.basis = m.group(1)

    # Basis functions
    m = _RE_BASIS_FUNCS.search(text)
    if m:
        result.nbasis = int(m.group(1))

    # Atom count from nuclear orientation table
    m = _RE_NATOMS_LINE.search(text)
    if m:
        result.natoms = len([l for l in m.group(1).strip().splitlines() if l.strip()])

    # SCF convergence
    scf_cycles = _RE_SCF_CYCLE.findall(text)
    if scf_cycles:
        result.scf.cycles = int(scf_cycles[-1][0])
    if _RE_SCF_CONVERGED.search(text):
        result.scf.converged = True

    # Energies
    m = _RE_SCF_ENERGY.search(text)
    if m:
        result.scf.energy = float(m.group(1))

    m = _RE_MP2_ENERGY.search(text)
    if m:
        result.mp2.energy = float(m.group(1))

    m = _RE_CCSD_CORR.search(text)
    if m:
        result.cc.ccsd_corr_energy = float(m.group(1))

    m = _RE_CCSD_TOTAL.search(text)
    if m:
        result.cc.ccsd_energy = float(m.group(1))
        result.cc.converged = True

    m = _RE_CCSD_T_TOTAL.search(text)
    if m:
        result.cc.ccsd_t_energy = float(m.group(1))

    m = _RE_CCSD_T_CORR.search(text)
    if m:
        result.cc.ccsd_t_corr_energy = float(m.group(1))

    if _RE_CC_CONVERGED.search(text):
        result.cc.converged = True

    # T diagnostic
    m = _RE_T_DIAG.search(text)
    if m:
        result.cc.t1_squared = float(m.group(1))
        result.cc.t2_squared = float(m.group(2))

    # Leading amplitudes
    for m in _RE_AMPLITUDE_DOUBLE.finditer(text):
        result.leading_amplitudes.append(Amplitude(
            value=float(m.group(1)),
            excitation_type="double",
            from_orbitals=m.group(2).strip(),
            to_orbitals=m.group(3).strip(),
        ))
    for m in _RE_AMPLITUDE_SINGLE.finditer(text):
        already = any(
            abs(a.value - float(m.group(1))) < 1e-10
            for a in result.leading_amplitudes
        )
        if not already:
            result.leading_amplitudes.append(Amplitude(
                value=float(m.group(1)),
                excitation_type="single",
                from_orbitals=m.group(2).strip(),
                to_orbitals=m.group(3).strip(),
            ))

    # Wall time
    m = _RE_TOTAL_TIME.search(text)
    if m:
        result.wall_time = m.group(1)

    # Check for .fchk
    fchk_path = path.with_suffix(".fchk")
    result.fchk_exists = fchk_path.exists()

    # Overall success = SCF converged and (for CC methods) CC converged
    is_cc = result.method.lower() in ("ccsd", "ccsd(t)")
    if is_cc:
        result.success = result.scf.converged and result.cc.converged
    else:
        result.success = result.scf.converged and result.mp2.energy is not None

    return result


# ---------------------------------------------------------------------------
# Batch parsing & reporting
# ---------------------------------------------------------------------------

def parse_batch(output_dir: str | Path) -> list[ParsedOutput]:
    """Parse all .out files in a directory."""
    output_dir = Path(output_dir)
    outs = sorted(output_dir.glob("*.out"))
    if not outs:
        print(f"  No .out files found in {output_dir}")
        return []

    results = [parse_output(f) for f in outs]
    return results


def print_summary(results: list[ParsedOutput]) -> None:
    """Print a human-readable summary table of parsed results."""
    if not results:
        print("  No results to display.")
        return

    header = f"{'Molecule':<50} {'OK':>3} {'Method':>8} {'Basis':>10} {'E(HF)':>16} {'E(MP2)':>16} {'E(CCSD)':>16} {'E(CCSD(T))':>16} {'T1²':>7} {'Time':>10}"
    print(header)
    print("─" * len(header))

    for r in results:
        name = r.filename.replace(".out", "")
        if len(name) > 49:
            name = name[:46] + "..."
        ok = "✓" if r.success else "✗"
        hf = f"{r.scf.energy:.8f}" if r.scf.energy else ""
        mp2 = f"{r.mp2.energy:.8f}" if r.mp2.energy else ""
        ccsd = f"{r.cc.ccsd_energy:.8f}" if r.cc.ccsd_energy else ""
        ccsd_t = f"{r.cc.ccsd_t_energy:.8f}" if r.cc.ccsd_t_energy else ""
        t1sq = f"{r.cc.t1_squared:.4f}" if r.cc.t1_squared is not None else ""
        wt = r.wall_time or ""
        print(f"{name:<50} {ok:>3} {r.method:>8} {r.basis:>10} {hf:>16} {mp2:>16} {ccsd:>16} {ccsd_t:>16} {t1sq:>7} {wt:>10}")

    succeeded = sum(1 for r in results if r.success)
    failed = [r for r in results if not r.success]
    print(f"\n  {succeeded}/{len(results)} succeeded")
    if failed:
        print("  Failed:")
        for r in failed:
            print(f"    {r.filename}: {r.error or 'SCF/CC did not converge'}")


def export_json(results: list[ParsedOutput], path: str | Path) -> None:
    """Export results as JSON."""
    path = Path(path)

    def _serialize(obj):
        if isinstance(obj, Path):
            return str(obj)
        raise TypeError(f"Cannot serialize {type(obj)}")

    data = []
    for r in results:
        d = asdict(r)
        data.append(d)

    path.write_text(json.dumps(data, indent=2, default=_serialize))
    print(f"  Results exported to {path}")


def export_csv(results: list[ParsedOutput], path: str | Path) -> None:
    """Export a flat energy summary as CSV."""
    path = Path(path)
    fieldnames = [
        "molecule", "success", "method", "basis", "natoms", "nbasis",
        "scf_energy", "mp2_energy", "ccsd_energy", "ccsd_t_energy",
        "ccsd_corr_energy", "t1_squared", "t2_squared",
        "scf_converged", "cc_converged", "fchk_exists", "wall_time",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({
                "molecule": r.filename.replace(".out", ""),
                "success": r.success,
                "method": r.method,
                "basis": r.basis,
                "natoms": r.natoms,
                "nbasis": r.nbasis,
                "scf_energy": r.scf.energy,
                "mp2_energy": r.mp2.energy,
                "ccsd_energy": r.cc.ccsd_energy,
                "ccsd_t_energy": r.cc.ccsd_t_energy,
                "ccsd_corr_energy": r.cc.ccsd_corr_energy,
                "t1_squared": r.cc.t1_squared,
                "t2_squared": r.cc.t2_squared,
                "scf_converged": r.scf.converged,
                "cc_converged": r.cc.converged,
                "fchk_exists": r.fchk_exists,
                "wall_time": r.wall_time,
            })
    print(f"  CSV exported to {path}")
