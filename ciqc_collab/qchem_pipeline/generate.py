"""Generate Q-Chem input files from XYZ coordinate files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import LARGE_ATOM_THRESHOLD, ProjectConfig, QChemParams


# ---------------------------------------------------------------------------
# Q-Chem input templates
# ---------------------------------------------------------------------------

def _build_rem_block(p: QChemParams, *, force_method: str | None = None) -> str:
    """Build the $rem section from parameters."""
    method = force_method or p.method
    is_mp2 = method.lower() == "mp2"

    lines = [
        f"    METHOD          {method}",
        f"    BASIS           {p.basis}",
    ]
    if p.mem_total:
        lines.append(f"    MEM_TOTAL       {p.mem_total}")
    lines += [
        f"    UNRESTRICTED    {'TRUE' if p.unrestricted else 'FALSE'}",
        f"    INTERNAL_STABILITY {'TRUE' if p.internal_stability else 'FALSE'}",
        f"    SCF_ALGORITHM   {p.scf_algorithm}",
        f"    MAX_SCF_CYCLES  {p.max_scf_cycles}",
    ]
    if not is_mp2:
        lines += [
            f"    THRESH          {p.thresh}",
            f"    SCF_CONVERGENCE {p.scf_convergence}",
            f"    CC_CONVERGENCE  {p.cc_convergence}",
        ]
    lines.append(f"    N_FROZEN_CORE   {p.n_frozen_core}")
    if p.gui:
        lines.append(f"    GUI             {p.gui}")
    return "\n".join(lines)


def build_input(
    coords: str,
    charge: int,
    mult: int,
    params: QChemParams,
    *,
    force_method: str | None = None,
) -> str:
    """Return a complete Q-Chem input file as a string."""
    rem = _build_rem_block(params, force_method=force_method)
    return (
        f"$molecule\n"
        f"{charge} {mult}\n"
        f"{coords}\n"
        f"$end\n"
        f"\n"
        f"$rem\n"
        f"{rem}\n"
        f"$end\n"
    )


# ---------------------------------------------------------------------------
# XYZ parsing
# ---------------------------------------------------------------------------

@dataclass
class MoleculeXYZ:
    name: str
    natoms: int
    comment: str
    coords: str
    source_path: Path


def read_xyz(path: str | Path) -> MoleculeXYZ:
    """Parse a standard XYZ file, returning atom coordinates (lines 3+)."""
    path = Path(path)
    lines = path.read_text().splitlines()
    natoms = int(lines[0].strip())
    comment = lines[1].strip() if len(lines) > 1 else ""
    coord_lines = lines[2 : 2 + natoms]
    coords = "\n".join(coord_lines)
    return MoleculeXYZ(
        name=path.stem,
        natoms=natoms,
        comment=comment,
        coords=coords,
        source_path=path.resolve(),
    )


# ---------------------------------------------------------------------------
# Batch generation
# ---------------------------------------------------------------------------

def _safe_basis_tag(basis: str) -> str:
    return basis.replace("(", "").replace(")", "").replace("*", "s").replace("+", "p").replace("/", "_")


def generate_inputs(
    xyz_dir: str | Path,
    output_dir: str | Path,
    config: ProjectConfig,
    *,
    method_override: str | None = None,
    basis_override: str | None = None,
    mem_override: int | None = None,
) -> list[Path]:
    """Generate Q-Chem .in files for every XYZ in *xyz_dir*.

    Returns the list of generated input file paths.
    """
    xyz_dir = Path(xyz_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    params = config.qchem
    if basis_override:
        params = QChemParams(**{**params.__dict__, "basis": basis_override})
    if mem_override is not None:
        params = QChemParams(**{**params.__dict__, "mem_total": mem_override})

    effective_method = method_override or params.method
    is_auto = effective_method.lower() == "auto"
    basis_tag = _safe_basis_tag(params.basis)

    xyz_files = sorted(xyz_dir.glob("*.xyz"))
    if not xyz_files:
        print(f"  No .xyz files found in {xyz_dir}")
        return []

    if is_auto:
        print(f"  Method: auto (≤{LARGE_ATOM_THRESHOLD} atoms → ccsd(t), >{LARGE_ATOM_THRESHOLD} atoms → mp2)\n")

    generated: list[Path] = []
    for xf in xyz_files:
        mol = read_xyz(xf)
        charge, mult = config.resolve_charge_mult(mol.name)
        is_large = mol.natoms > LARGE_ATOM_THRESHOLD

        if is_auto:
            method_for_mol = "mp2" if is_large else "ccsd(t)"
        elif is_large and effective_method.lower() not in ("mp2",):
            method_for_mol = "mp2"
            print(f"  {mol.name}: {mol.natoms} atoms — too large for {effective_method}, forcing MP2")
        else:
            method_for_mol = effective_method

        method_tag = method_for_mol.replace("(", "").replace(")", "")
        tag = f"{method_tag}_{basis_tag}"
        out_name = f"{mol.name}_{tag}.in"
        content = build_input(mol.coords, charge, mult, params, force_method=method_for_mol)

        out_path = output_dir / out_name
        out_path.write_text(content)
        generated.append(out_path)
        print(f"  {out_name}  ({mol.natoms} atoms, {method_for_mol}, charge={charge}, mult={mult})")

    print(f"\n  Generated {len(generated)} input files in {output_dir}")
    return generated
