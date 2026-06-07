from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch.utils.data import Dataset


def _parse_fchk_mo_data(fchk_path: Path) -> dict[str, np.ndarray]:
    """Parse .fchk for MO coefficients and AO metadata.

    Returns dict with:
        mo_coeffs: (norb, nbasis) -- MO coefficient matrix (columns = MOs)
        ao_nuclear_charges: (nbasis,) -- nuclear charge Z of each AO's atom
        ao_angular_momentum: (nbasis,) -- angular momentum l of each AO
    """
    lines = fchk_path.read_text().splitlines()

    def _read_array(start: int, n: int, dtype=float):
        vals: list = []
        i = start + 1
        while len(vals) < n:
            vals.extend(dtype(x) for x in lines[i].split())
            i += 1
        return vals[:n]

    nbasis = 0
    atomic_nums: list[int] = []
    shell_types: list[int] = []
    shell_to_atom: list[int] = []
    mo_coeffs_flat: list[float] = []

    for i, line in enumerate(lines):
        if "Number of basis functions" in line:
            nbasis = int(line.split("I")[1])
        elif "Atomic numbers" in line and "N=" in line:
            n = int(line.split("N=")[1])
            atomic_nums = _read_array(i, n, int)
        elif line.startswith("Shell types") and "N=" in line:
            n = int(line.split("N=")[1])
            shell_types = _read_array(i, n, int)
        elif line.startswith("Shell to atom map") and "N=" in line:
            n = int(line.split("N=")[1])
            shell_to_atom = _read_array(i, n, int)
        elif "Alpha MO coefficients" in line:
            n = int(line.split("N=")[1])
            mo_coeffs_flat = _read_array(i, n, float)

    ao_Z: list[int] = []
    ao_l: list[int] = []
    for stype, atom_idx in zip(shell_types, shell_to_atom):
        Z = atomic_nums[atom_idx - 1]
        if stype == -1:  # SP shell
            ao_Z.append(Z); ao_l.append(0)
            for _ in range(3):
                ao_Z.append(Z); ao_l.append(1)
        elif stype == 0:
            ao_Z.append(Z); ao_l.append(0)
        elif stype == 1:
            for _ in range(3):
                ao_Z.append(Z); ao_l.append(1)
        else:
            n_funcs = 2 * abs(stype) + 1 if stype > 0 else (abs(stype) * 2 + 1 + 1)
            if stype < -1:
                n_funcs = (abs(stype) + 1) * (abs(stype) + 2) // 2
            for _ in range(n_funcs):
                ao_Z.append(Z); ao_l.append(abs(stype))

    C = np.array(mo_coeffs_flat).reshape(nbasis, nbasis)

    return {
        "mo_coeffs": C.T.astype(np.float32),
        "ao_nuclear_charges": np.array(ao_Z, dtype=np.float32),
        "ao_angular_momentum": np.array(ao_l, dtype=np.float32),
    }


def _parse_orbital_energies(path: Path) -> np.ndarray:
    """Parse mo_ene_for_qis.dat: occupied then virtual orbital energies."""
    with open(path) as f:
        lines = f.readlines()

    occ_energies: list[float] = []
    virt_energies: list[float] = []
    section = None

    for line in lines:
        line = line.strip()
        if line == "O":
            section = "occ"
            continue
        elif line == "V":
            section = "virt"
            continue
        if section is None:
            continue
        tokens = line.split()
        if len(tokens) == 2 and tokens[0] == "1":
            continue
        target = occ_energies if section == "occ" else virt_energies
        target.extend(float(x) for x in tokens)

    return np.concatenate([
        np.array(occ_energies, dtype=np.float64),
        np.array(virt_energies, dtype=np.float64),
    ])


def _parse_amplitude_file(path: Path) -> np.ndarray:
    with open(path) as f:
        lines = f.readlines()

    shape_info = lines[1].split()
    ndim = int(shape_info[0])
    shape = tuple(int(x) for x in shape_info[1 : ndim + 1])

    values: list[float] = []
    for line in lines[2:]:
        values.extend(float(x) for x in line.split())

    arr = np.array(values, dtype=np.float64)
    expected = int(np.prod(shape))
    if arr.size != expected:
        raise ValueError(
            f"Expected {expected} values for shape {shape}, got {arr.size} in {path}"
        )
    return arr.reshape(shape)


class CCSDAmplitudeDataset(Dataset):
    """Load precomputed CCSD amplitudes (inputs) and LUCJ parameters (targets).

    Each valid job directory must contain both `.done` and `.lucj_done` marker
    files. The list of valid jobs is cached at init; tensors are loaded lazily
    in __getitem__.
    """

    def __init__(self, jobs_dir: str | Path, species_filter: str | None = None) -> None:
        """
        Args:
            jobs_dir: Path to directory containing job subdirectories.
            species_filter: Filter jobs by species type. Options:
                None  — all species (R, TS, P)
                "TS"  — only transition states
                "RP"  — only reactants and products (no TS)
        """
        self.jobs_dir = Path(jobs_dir)
        self.job_names: list[str] = []
        self._norbs: list[int] = []

        _allowed_suffixes: set[str] | None = None
        if species_filter == "TS":
            _allowed_suffixes = {"_TS"}
        elif species_filter == "RP":
            _allowed_suffixes = {"_R", "_P"}

        for d in sorted(self.jobs_dir.iterdir()):
            if not d.is_dir():
                continue
            if (d / ".done").exists() and (d / ".lucj_done").exists():
                if _allowed_suffixes is not None:
                    suffix = "_" + d.name.rsplit("_", 1)[-1]
                    if suffix not in _allowed_suffixes:
                        continue
                self.job_names.append(d.name)
                meta_path = d / "lucj_metadata.json"
                with open(meta_path) as f:
                    self._norbs.append(json.load(f)["norb"])

    def get_norb(self, idx: int) -> int:
        return self._norbs[idx]

    def __len__(self) -> int:
        return len(self.job_names)

    def __getitem__(self, idx: int) -> dict[str, Any]:
        name = self.job_names[idx]
        job_dir = self.jobs_dir / name

        t1 = _parse_amplitude_file(job_dir / "ccsd_t1.dat")
        t2 = _parse_amplitude_file(job_dir / "ccsd_t2.dat")

        with open(job_dir / "lucj_metadata.json") as f:
            meta = json.load(f)

        norb: int = meta["norb"]
        n_reps: int = meta["n_reps"]
        nocc: int = meta["nocc"]
        nvirt: int = meta["nvirt"]

        J = np.load(job_dir / "lucj_diag_coulomb_mats.npy")
        U = np.load(job_dir / "lucj_orbital_rotations.npy")
        final_rot = np.load(job_dir / "lucj_final_orbital_rotation.npy")

        kappa_real = np.load(job_dir / "kappa_real.npy")
        kappa_imag = np.load(job_dir / "kappa_imag.npy")

        orb_ene_path = job_dir / "mo_ene_for_qis.dat"
        if orb_ene_path.exists():
            orb_energies = _parse_orbital_energies(orb_ene_path).astype(np.float32)
        else:
            orb_energies = np.zeros(norb, dtype=np.float32)

        fchk_files = list(job_dir.glob("*.fchk"))
        if fchk_files:
            mo_data = _parse_fchk_mo_data(fchk_files[0])
            nbasis_fchk = mo_data["mo_coeffs"].shape[0]
            if nbasis_fchk < norb:
                mo_coeffs_padded = np.zeros((norb, mo_data["mo_coeffs"].shape[1]),
                                            dtype=np.float32)
                mo_coeffs_padded[:nbasis_fchk] = mo_data["mo_coeffs"]
                mo_coeffs_padded[nbasis_fchk:norb] = mo_data["mo_coeffs"][:norb - nbasis_fchk]
                mo_data["mo_coeffs"] = mo_coeffs_padded
        else:
            nbasis = norb
            mo_data = {
                "mo_coeffs": np.zeros((norb, nbasis), dtype=np.float32),
                "ao_nuclear_charges": np.zeros(nbasis, dtype=np.float32),
                "ao_angular_momentum": np.zeros(nbasis, dtype=np.float32),
            }

        return {
            "t1": torch.from_numpy(t1.astype(np.float32)),
            "t2": torch.from_numpy(t2.astype(np.float32)),
            "J": torch.from_numpy(J.astype(np.float32)),
            "U_real": torch.from_numpy(U.real.astype(np.float32)),
            "U_imag": torch.from_numpy(U.imag.astype(np.float32)),
            "kappa_real": torch.from_numpy(kappa_real.astype(np.float32)),
            "kappa_imag": torch.from_numpy(kappa_imag.astype(np.float32)),
            "final_rot": torch.from_numpy(final_rot.astype(np.float32)),
            "orb_energies": torch.from_numpy(orb_energies),
            "mo_coeffs": torch.from_numpy(mo_data["mo_coeffs"]),
            "ao_Z": torch.from_numpy(mo_data["ao_nuclear_charges"]),
            "ao_l": torch.from_numpy(mo_data["ao_angular_momentum"]),
            "nocc": nocc,
            "nvirt": nvirt,
            "norb": norb,
            "nbasis": mo_data["mo_coeffs"].shape[1],
            "n_reps": n_reps,
            "name": name,
        }

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(jobs_dir={self.jobs_dir!r}, n_jobs={len(self)})"

    @staticmethod
    def collate_fn(batch: list[dict[str, Any]]) -> dict[str, Any]:
        max_nocc = max(s["nocc"] for s in batch)
        max_nvirt = max(s["nvirt"] for s in batch)
        max_norb = max(s["norb"] for s in batch)
        max_n_reps = max(s["n_reps"] for s in batch)
        max_nbasis = max(s["nbasis"] for s in batch)
        B = len(batch)

        t1 = torch.zeros(B, max_nocc, max_nvirt)
        t2 = torch.zeros(B, max_nocc, max_nocc, max_nvirt, max_nvirt)
        J = torch.zeros(B, max_n_reps, 2, max_norb, max_norb)
        U_real = torch.zeros(B, max_n_reps, max_norb, max_norb)
        U_imag = torch.zeros(B, max_n_reps, max_norb, max_norb)
        kappa_real = torch.zeros(B, max_n_reps, max_norb, max_norb)
        kappa_imag = torch.zeros(B, max_n_reps, max_norb, max_norb)
        final_rot = torch.zeros(B, max_norb, max_norb)
        orb_energies = torch.zeros(B, max_norb)
        mo_coeffs = torch.zeros(B, max_norb, max_nbasis)
        ao_Z = torch.zeros(B, max_nbasis)
        ao_l = torch.zeros(B, max_nbasis)

        orb_mask = torch.zeros(B, max_norb, dtype=torch.bool)
        occ_mask = torch.zeros(B, max_nocc, dtype=torch.bool)
        virt_mask = torch.zeros(B, max_nvirt, dtype=torch.bool)

        nocc_list: list[int] = []
        nvirt_list: list[int] = []
        norb_list: list[int] = []
        n_reps_list: list[int] = []
        names: list[str] = []

        for i, s in enumerate(batch):
            no, nv, nb, nr = s["nocc"], s["nvirt"], s["norb"], s["n_reps"]

            nbas = s["nbasis"]

            t1[i, :no, :nv] = s["t1"]
            t2[i, :no, :no, :nv, :nv] = s["t2"]
            J[i, :nr, :, :nb, :nb] = s["J"]
            U_real[i, :nr, :nb, :nb] = s["U_real"]
            U_imag[i, :nr, :nb, :nb] = s["U_imag"]
            kappa_real[i, :nr, :nb, :nb] = s["kappa_real"]
            kappa_imag[i, :nr, :nb, :nb] = s["kappa_imag"]
            final_rot[i, :nb, :nb] = s["final_rot"]
            orb_energies[i, :nb] = s["orb_energies"]
            mo_coeffs[i, :nb, :nbas] = s["mo_coeffs"]
            ao_Z[i, :nbas] = s["ao_Z"]
            ao_l[i, :nbas] = s["ao_l"]

            orb_mask[i, :nb] = True
            occ_mask[i, :no] = True
            virt_mask[i, :nv] = True

            nocc_list.append(no)
            nvirt_list.append(nv)
            norb_list.append(nb)
            n_reps_list.append(nr)
            names.append(s["name"])

        return {
            "t1": t1,
            "t2": t2,
            "orb_energies": orb_energies,
            "mo_coeffs": mo_coeffs,
            "ao_Z": ao_Z,
            "ao_l": ao_l,
            "j_target": J,
            "u_real_target": U_real,
            "u_imag_target": U_imag,
            "kappa_real_target": kappa_real,
            "kappa_imag_target": kappa_imag,
            "final_rot": final_rot,
            "orb_mask": orb_mask,
            "occ_mask": occ_mask,
            "virt_mask": virt_mask,
            "nocc": nocc_list,
            "nvirt": nvirt_list,
            "norbs": torch.tensor(norb_list),
            "n_reps": n_reps_list,
            "names": names,
        }
