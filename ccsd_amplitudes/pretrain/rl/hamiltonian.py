"""Per-molecule active-space Hamiltonians for the energy reward.

Rebuilds the exact pyscf RHF used by generate_rhf_dataset.py (STO-3G, frozen
core = number of non-H atoms) and stores the ffsim MolecularHamiltonian
tensors in the active space so that reward workers never touch pyscf again:

    <ham-dir>/<name>.npz : one_body (n,n), two_body (n,n,n,n), constant,
                           norb, nelec_a, nelec_b, e_hf, e_ccsd, nfrozen

For norb=44 the two-body tensor is 44^4*8 B = 30 MB; norb=29 -> 5.7 MB.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

_Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
      "F": 9, "Ne": 10, "P": 15, "S": 16, "Cl": 17}


def read_xyz(path: Path) -> str:
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    return "\n".join(lines[2:2 + n])


def n_frozen_core(atom_block: str) -> int:
    return sum(1 for ln in atom_block.splitlines() if _Z.get(ln.split()[0], 0) >= 3)


def build_one(args) -> str:
    name, jobs_dir, ham_dir = args
    out = Path(ham_dir) / f"{name}.npz"
    if out.exists():
        return "cached"
    try:
        import ffsim
        from pyscf import cc, gto, scf

        xyz = read_xyz(Path(jobs_dir) / name / f"{name}.xyz")
        mol = gto.M(atom=xyz, basis="sto-3g", verbose=0)
        mf = scf.RHF(mol).run()
        if not mf.converged:
            return "hf_unconverged"
        nfrozen = n_frozen_core(xyz)
        active = list(range(nfrozen, mol.nao))
        mol_data = ffsim.MolecularData.from_scf(mf, active_space=active)
        ham = mol_data.hamiltonian
        mycc = cc.CCSD(mf, frozen=nfrozen).run()
        np.savez(
            out,
            one_body=ham.one_body_tensor.real.astype(np.float64),
            two_body=ham.two_body_tensor.real.astype(np.float64),
            constant=float(ham.constant),
            norb=int(mol_data.norb), nelec_a=int(mol_data.nelec[0]),
            nelec_b=int(mol_data.nelec[1]), nfrozen=nfrozen,
            e_hf=float(mf.e_tot), e_ccsd=float(mycc.e_tot),
            ccsd_converged=bool(mycc.converged),
        )
        return "ok"
    except Exception as e:  # noqa: BLE001
        return f"err:{type(e).__name__}:{e}"


def load_hamiltonian(ham_dir: str | Path, name: str):
    """Returns (ffsim.MolecularHamiltonian, norb, nelec, e_hf, e_ccsd)."""
    import ffsim
    d = np.load(Path(ham_dir) / f"{name}.npz")
    ham = ffsim.MolecularHamiltonian(
        one_body_tensor=d["one_body"], two_body_tensor=d["two_body"],
        constant=float(d["constant"]))
    return (ham, int(d["norb"]), (int(d["nelec_a"]), int(d["nelec_b"])),
            float(d["e_hf"]), float(d["e_ccsd"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names-file", required=True)
    ap.add_argument("--jobs-dir", default="jobs")
    ap.add_argument("--ham-dir", default="rhf_hamiltonians")
    ap.add_argument("--n-procs", type=int, default=16)
    args = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    Path(args.ham_dir).mkdir(exist_ok=True)
    names = [ln.strip() for ln in open(args.names_file) if ln.strip()]
    t0 = time.time()
    counts: dict[str, int] = {}
    with Pool(args.n_procs) as pool:
        for i, r in enumerate(pool.imap_unordered(
                build_one, [(n, args.jobs_dir, args.ham_dir) for n in names])):
            key = r.split(":")[0]
            counts[key] = counts.get(key, 0) + 1
            if r.startswith("err") and counts[key] <= 5:
                print("  ", r, flush=True)
            if (i + 1) % 100 == 0:
                print(f"  {i+1}/{len(names)}  {counts}  {time.time()-t0:.0f}s", flush=True)
    print(json.dumps(counts), f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
