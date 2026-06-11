#!/usr/bin/env python3
"""Generate the RHF-CCSD dataset for the reconstruction-loss surrogate.

For each molecule (geometry read from jobs/<name>/<name>.xyz), runs
pyscf RHF-CCSD/STO-3G with frozen core and stores the restricted spatial
amplitudes + orbital info needed by the model:

  t1            (nocc, nvirt)
  t2            (nocc, nocc, nvirt, nvirt)   <- ffsim's expected restricted t2
  orb_energies  (norb,)                       active-space RHF eigenvalues
  mo_coeff      (nao, norb)                    active-space MO coefficients
  ao_Z          (nao,)                         nuclear charge of each AO's atom
  ao_l          (nao,)                         angular momentum of each AO
  nocc, nvirt, norb, nfrozen, e_hf, e_ccsd

Pure reconstruction training only needs t2; the rest supports the orbital
encodings and diagnostics. No ffsim here (the model outputs its own U,Z).

Shardable for parallelism:  --shard i --n-shards N

Usage:
    python generate_rhf_dataset.py --shard 0 --n-shards 16
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

_L = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
_Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
      "F": 9, "Ne": 10, "P": 15, "S": 16, "Cl": 17}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    return "\n".join(lines[2:2 + n])


def n_frozen_core(atom_block):
    return sum(1 for ln in atom_block.splitlines()
               if _Z.get(ln.split()[0], 0) >= 3)


def ao_info(mol):
    """Per-AO nuclear charge and angular momentum from ao_labels."""
    aoZ, aoL = [], []
    charges = mol.atom_charges()
    for lbl in mol.ao_labels(fmt=False):
        atom_idx, _sym, nl, _m = lbl
        aoZ.append(int(charges[atom_idx]))
        aoL.append(_L.get(nl[-1], 0))
    return np.array(aoZ, np.int16), np.array(aoL, np.int16)


def process(name, jobs_dir, out_dir):
    from pyscf import gto, scf, cc
    xyz_path = jobs_dir / name / f"{name}.xyz"
    if not xyz_path.exists():
        return "no_xyz"
    out_path = out_dir / f"{name}.npz"
    if out_path.exists():
        return "cached"

    xyz = read_xyz(xyz_path)
    mol = gto.M(atom=xyz, basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    if not mf.converged:
        return "hf_unconverged"
    nfrozen = n_frozen_core(xyz)
    mycc = cc.CCSD(mf, frozen=nfrozen)
    mycc.kernel()
    if not mycc.converged:
        return "ccsd_unconverged"

    t1, t2 = mycc.t1, mycc.t2
    nocc, nvirt = t1.shape
    norb = nocc + nvirt
    aoZ, aoL = ao_info(mol)

    np.savez(
        out_path,
        t1=t1.astype(np.float32),
        t2=t2.astype(np.float32),
        orb_energies=mf.mo_energy[nfrozen:].astype(np.float32),
        mo_coeff=mf.mo_coeff[:, nfrozen:].astype(np.float32),
        ao_Z=aoZ, ao_l=aoL,
        nocc=nocc, nvirt=nvirt, norb=norb, nfrozen=nfrozen,
        e_hf=float(mf.e_tot), e_ccsd=float(mycc.e_tot),
    )
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs-dir", default="jobs")
    ap.add_argument("--out-dir", default="rhf_dataset")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    args = ap.parse_args()

    jobs_dir = Path(args.jobs_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)
    manifest = json.load(open(jobs_dir / "manifest.json"))

    names = []
    for e in manifest:
        for sp in ("R", "TS", "P"):
            names.append(f"{e['formula']}_{e['rxn_id']}_{sp}")
    names = names[args.shard::args.n_shards]

    counts = {}
    t0 = time.time()
    for i, name in enumerate(names):
        try:
            r = process(name, jobs_dir, out_dir)
        except Exception as e:
            r = f"err:{type(e).__name__}"
        counts[r] = counts.get(r, 0) + 1
        if (i + 1) % 50 == 0:
            done = sum(v for k, v in counts.items() if k in ("ok", "cached"))
            print(f"[shard {args.shard}] {i+1}/{len(names)}  ok+cached={done}  "
                  f"{counts}  {time.time()-t0:.0f}s", flush=True)

    print(f"[shard {args.shard}] DONE {counts}  {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
