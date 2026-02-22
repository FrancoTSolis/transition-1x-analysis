"""Molecule registry, Q-Chem defaults, and YAML configuration loading."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


# ---------------------------------------------------------------------------
# Molecule metadata
# ---------------------------------------------------------------------------

KNOWN_CHARGE_MULT: dict[str, tuple[int, int]] = {
    "02_Cationic_Aza_Cope": (1, 1),
}
DEFAULT_CHARGE_MULT = (0, 1)

LARGE_ATOM_THRESHOLD = 40


def charge_mult_for(filename: str) -> tuple[int, int]:
    """Return (charge, multiplicity) for a molecule by matching filename prefix."""
    for prefix, cm in KNOWN_CHARGE_MULT.items():
        if filename.startswith(prefix):
            return cm
    return DEFAULT_CHARGE_MULT


# ---------------------------------------------------------------------------
# Q-Chem defaults (matching collaborator recommendations)
# ---------------------------------------------------------------------------

QCHEM_DEFAULTS: dict[str, Any] = {
    "unrestricted": True,
    "internal_stability": True,
    "scf_algorithm": "DIIS_GDM",
    "max_scf_cycles": 200,
    "thresh": 15,
    "scf_convergence": 10,
    "cc_convergence": 7,
    "n_frozen_core": "FC",
    "gui": 2,
}

SLURM_DEFAULTS: dict[str, Any] = {
    "account": None,
    "time": "08:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 18,
    "constraint": "cpu",
    "qos": "regular",
    "modules": ["qchem"],
}


# ---------------------------------------------------------------------------
# Project configuration (loaded from YAML)
# ---------------------------------------------------------------------------

@dataclass
class QChemParams:
    method: str = "auto"
    basis: str = "STO-3G"
    mem_total: int | None = None
    unrestricted: bool = True
    internal_stability: bool = True
    scf_algorithm: str = "DIIS_GDM"
    max_scf_cycles: int = 200
    thresh: int = 15
    scf_convergence: int = 10
    cc_convergence: int = 7
    n_frozen_core: str = "FC"
    gui: int = 2


@dataclass
class SlurmParams:
    account: str | None = None
    time: str = "08:00:00"
    nodes: int = 1
    ntasks: int = 1
    cpus_per_task: int = 18
    constraint: str = "cpu"
    qos: str = "regular"
    modules: list[str] = field(default_factory=lambda: ["qchem"])
    extra_exports: dict[str, str] = field(default_factory=dict)


@dataclass
class ProjectConfig:
    qchem: QChemParams = field(default_factory=QChemParams)
    slurm: SlurmParams = field(default_factory=SlurmParams)
    charge_mult_overrides: dict[str, tuple[int, int]] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ProjectConfig":
        path = Path(path)
        if not path.exists():
            return cls()
        with open(path) as f:
            raw = yaml.safe_load(f) or {}

        qc_raw = raw.get("qchem", {})
        slurm_raw = raw.get("slurm", {})
        cm_raw = raw.get("charge_mult_overrides", {})

        qchem = QChemParams(**{k: v for k, v in qc_raw.items() if k in QChemParams.__dataclass_fields__})
        slurm = SlurmParams(**{k: v for k, v in slurm_raw.items() if k in SlurmParams.__dataclass_fields__})
        overrides = {k: tuple(v) for k, v in cm_raw.items()}

        return cls(qchem=qchem, slurm=slurm, charge_mult_overrides=overrides)

    def resolve_charge_mult(self, filename: str) -> tuple[int, int]:
        for prefix, cm in self.charge_mult_overrides.items():
            if filename.startswith(prefix):
                return cm
        return charge_mult_for(filename)


def find_config(start_dir: str | Path | None = None) -> Path | None:
    """Walk up from *start_dir* looking for ``pipeline_config.yaml``."""
    d = Path(start_dir or os.getcwd()).resolve()
    for _ in range(10):
        candidate = d / "pipeline_config.yaml"
        if candidate.exists():
            return candidate
        if d.parent == d:
            break
        d = d.parent
    return None
