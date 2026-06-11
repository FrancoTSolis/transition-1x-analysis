#!/usr/bin/env python3
"""Build a cache of RHF-CCSD + LUCJ results for learnability analysis.

Runs pyscf RHF-CCSD/STO-3G (frozen core) on molecules, computes LUCJ via
ffsim, and saves t1/t2/kappa/J per molecule to rhf_cache/{name}.npz.
Incremental: skips already-cached molecules.

Usage:
    python -m pretrain.build_rhf_cache --target-norb 22 --n-mols 120
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    natoms = int(lines[0])
    return "\n".join(lines[2 : 2 + natoms])


def n_frozen_core(atom_block):
    frozen = 0
    for line in atom_block.splitlines():
        sym = line.split()[0]
        Z = {"H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "F": 9,
             "S": 16, "Cl": 17, "P": 15}.get(sym, 0)
        if Z >= 3:
            frozen += 1
    return frozen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-norb", type=int, default=22)
    parser.add_argument("--n-mols", type=int, default=120)
    parser.add_argument("--jobs-dir", default="jobs")
    parser.add_argument("--cache-dir", default="rhf_cache")
    parser.add_argument("--n-reps", type=int, default=2)
    parser.add_argument("--maxiter", type=int, default=30)
    parser.add_argument("--t2-only", action="store_true",
                        help="Only compute/cache RHF-CCSD t2/t1 (+orbitals); skip the "
                             "expensive LUCJ. Fast: ~5s/molecule.")
    args = parser.parse_args()

    from pyscf import gto, scf, cc
    import ffsim
    from scipy.linalg import logm

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(exist_ok=True)
    jobs_dir = Path(args.jobs_dir)
    manifest = json.load(open(jobs_dir / "manifest.json"))

    names = []
    for entry in manifest:
        for sp in ("R", "TS", "P"):
            names.append(f"{entry['formula']}_{entry['rxn_id']}_{sp}")

    n_cached = 0
    t_start = time.time()

    for name in names:
        if n_cached >= args.n_mols:
            break
        out_path = cache_dir / f"{name}.npz"
        if out_path.exists():
            n_cached += 1
            continue

        xyz_path = jobs_dir / name / f"{name}.xyz"
        if not xyz_path.exists():
            continue
        xyz = read_xyz(xyz_path)

        try:
            mol = gto.M(atom=xyz, basis="sto-3g", verbose=0)
            mf = scf.RHF(mol).run()
            nfrozen = n_frozen_core(xyz)
            mycc = cc.CCSD(mf, frozen=nfrozen)
            mycc.kernel()
            if not mycc.converged:
                continue
            t1, t2 = mycc.t1, mycc.t2
            nocc, nvirt = t1.shape
            norb = nocc + nvirt
            if norb != args.target_norb:
                continue

            orb_energies = mf.mo_energy[nfrozen:]
            if args.t2_only:
                np.savez(
                    out_path,
                    t1=t1.astype(np.float32),
                    t2=t2.astype(np.float32),
                    orb_energies=orb_energies.astype(np.float32),
                    nocc=nocc, nvirt=nvirt, norb=norb,
                )
                n_cached += 1
                if n_cached % 20 == 0:
                    print(f"  cached {n_cached} (t2-only, {time.time()-t_start:.0f}s)",
                          flush=True)
                continue

            pairs_aa = [(p, p + 1) for p in range(norb - 1)]
            pairs_ab = [(p, p) for p in range(norb)]
            op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
                t2=t2, t1=t1, n_reps=args.n_reps,
                interaction_pairs=(pairs_aa, pairs_ab),
                optimize=True, options=dict(maxiter=args.maxiter),
            )
            kappa0 = logm(op.orbital_rotations[0])
            kappa1 = logm(op.orbital_rotations[1])
            orb_energies = mf.mo_energy[nfrozen:]  # active orbital energies

            np.savez(
                out_path,
                t1=t1.astype(np.float32),
                t2=t2.astype(np.float32),
                kappa_real=np.stack([kappa0.real, kappa1.real]).astype(np.float32),
                kappa_imag=np.stack([kappa0.imag, kappa1.imag]).astype(np.float32),
                J=op.diag_coulomb_mats.astype(np.float32),
                # raw DF tensors for exact reconstruction:
                # orbital_rotations (U): (n_reps, norb, norb) complex
                # Z (single diag-coulomb per rep): (n_reps, norb, norb) real
                U_real=op.orbital_rotations.real.astype(np.float32),
                U_imag=op.orbital_rotations.imag.astype(np.float32),
                Z=op.diag_coulomb_mats[:, 0].astype(np.float32),
                orb_energies=orb_energies.astype(np.float32),
                nocc=nocc, nvirt=nvirt, norb=norb,
            )
            n_cached += 1
            if n_cached % 10 == 0:
                print(f"  cached {n_cached} (norb={norb}, {time.time()-t_start:.0f}s)",
                      flush=True)
        except Exception as e:
            print(f"  SKIP {name}: {e}", flush=True)
            continue

    print(f"\nCache built: {n_cached} molecules at norb={args.target_norb} "
          f"in {time.time()-t_start:.0f}s -> {cache_dir}/", flush=True)


if __name__ == "__main__":
    main()
