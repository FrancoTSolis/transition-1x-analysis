#!/usr/bin/env python3
"""Validate the RHF-CCSD fix vs the spin-orbital approach.

For a sample of molecules:
1. Run fresh pyscf RHF-CCSD/STO-3G (frozen core) -> restricted spatial t2
2. Compute LUCJ (kappa, J) from the correct RHF t2 via ffsim
3. Run a linear baseline: does RHF-t2 -> kappa/J have positive R²?

This tests whether the corrected (RHF restricted) pipeline produces a
learnable t2 -> LUCJ mapping, justifying (or not) a full recompute.

Usage:
    python -m pretrain.validate_rhf --n-mols 60 --formula-filter C2
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    natoms = int(lines[0])
    atom_lines = lines[2 : 2 + natoms]
    return "\n".join(atom_lines)


def n_frozen_core(atom_block):
    """Q-Chem FC: freeze 1s of each non-H first/second-row atom."""
    frozen = 0
    for line in atom_block.splitlines():
        sym = line.split()[0]
        Z = {"H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "F": 9,
             "S": 16, "Cl": 17, "P": 15}.get(sym, 0)
        if Z >= 3:  # Li and beyond freeze 1s
            frozen += 1
    return frozen


def run_rhf_ccsd(xyz_block, basis="sto-3g"):
    from pyscf import gto, scf, cc
    mol = gto.M(atom=xyz_block, basis=basis, verbose=0)
    mf = scf.RHF(mol).run()
    nfrozen = n_frozen_core(xyz_block)
    mycc = cc.CCSD(mf, frozen=nfrozen)
    mycc.kernel()
    return mycc.t1, mycc.t2, mycc.converged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-mols", type=int, default=60)
    parser.add_argument("--jobs-dir", default="jobs")
    parser.add_argument("--n-reps", type=int, default=2)
    parser.add_argument("--maxiter", type=int, default=50)
    args = parser.parse_args()

    import ffsim
    from scipy.linalg import logm
    from sklearn.linear_model import Ridge
    from sklearn.metrics import r2_score

    jobs_dir = Path(args.jobs_dir)
    manifest = json.load(open(jobs_dir.parent / "jobs" / "manifest.json"))

    # Group by spatial norb so we can run a fixed-size linear baseline
    samples = []  # (t1_rhf, t2_rhf, kappa_real, kappa_imag, J, norb, nocc, nvirt)
    t_start = time.time()

    names_to_try = []
    for entry in manifest:
        for sp in ("R", "TS", "P"):
            names_to_try.append(f"{entry['formula']}_{entry['rxn_id']}_{sp}")

    for name in names_to_try:
        if len(samples) >= args.n_mols:
            break
        xyz_path = jobs_dir / name / f"{name}.xyz"
        if not xyz_path.exists():
            continue
        xyz = read_xyz(xyz_path)
        try:
            t1, t2, conv = run_rhf_ccsd(xyz)
        except Exception as e:
            print(f"  SKIP {name}: {e}")
            continue
        if not conv:
            print(f"  SKIP {name}: CCSD not converged")
            continue

        nocc, nvirt = t1.shape
        norb = nocc + nvirt

        pairs_aa = [(p, p + 1) for p in range(norb - 1)]
        pairs_ab = [(p, p) for p in range(norb)]
        try:
            op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
                t2=t2, t1=t1, n_reps=args.n_reps,
                interaction_pairs=(pairs_aa, pairs_ab),
                optimize=True, options=dict(maxiter=args.maxiter),
            )
        except Exception as e:
            print(f"  SKIP {name}: ffsim {e}")
            continue

        kappa = logm(op.orbital_rotations[0])
        samples.append({
            "name": name, "t1": t1, "t2": t2,
            "kr": kappa.real, "ki": kappa.imag,
            "J": op.diag_coulomb_mats[0],
            "norb": norb, "nocc": nocc, "nvirt": nvirt,
        })
        if len(samples) % 10 == 0:
            print(f"  {len(samples)} done ({time.time()-t_start:.0f}s)")

    print(f"\nCollected {len(samples)} RHF-CCSD + LUCJ samples in {time.time()-t_start:.0f}s")

    # Group by norb
    from collections import defaultdict
    by_norb = defaultdict(list)
    for s in samples:
        by_norb[s["norb"]].append(s)

    print(f"norb distribution: {sorted((k, len(v)) for k, v in by_norb.items())}")

    # Pick the largest group for the linear baseline
    best_norb = max(by_norb, key=lambda k: len(by_norb[k]))
    group = by_norb[best_norb]
    if len(group) < 20:
        print(f"\nLargest group (norb={best_norb}) has only {len(group)} samples — "
              f"need >=20 for a meaningful baseline. Run with more --n-mols.")
        return

    nocc = group[0]["nocc"]
    nvirt = group[0]["nvirt"]
    norb = best_norb
    n = len(group)
    n_train = int(0.8 * n)
    print(f"\n=== Linear baseline on norb={norb} (nocc={nocc}, nvirt={nvirt}), {n} samples ===")

    # Features: t2 summary
    feats = np.array([
        np.concatenate([
            s["t2"].mean(axis=(2, 3)).ravel(),
            s["t2"].mean(axis=(0, 1)).ravel(),
            [s["t2"].std(), np.abs(s["t2"]).max()],
        ]) for s in group
    ])
    Xtr, Xte = feats[:n_train], feats[n_train:]

    triu = np.triu_indices(norb, k=1)
    results = []

    # J_aa off-diag
    y = np.array([[s["J"][0, i, i + 1] for i in range(norb - 1)] for s in group])
    r = Ridge(1.0).fit(Xtr, y[:n_train])
    results.append(("J_aa off-diag", r2_score(y[n_train:], r.predict(Xte))))

    # kappa_real upper tri
    y = np.array([s["kr"][triu] for s in group])
    r = Ridge(10.0).fit(Xtr, y[:n_train])
    results.append(("kappa_re upper tri", r2_score(y[n_train:].ravel(), r.predict(Xte).ravel())))

    # kappa_real occ-virt block
    y = np.array([s["kr"][:nocc, nocc:].ravel() for s in group])
    r = Ridge(10.0).fit(Xtr, y[:n_train])
    results.append(("kappa_re occ-virt", r2_score(y[n_train:].ravel(), r.predict(Xte).ravel())))

    print(f"\n{'Target':<25} {'R²':>8} {'Signal':>8}")
    print("-" * 45)
    for lbl, r2v in results:
        sig = "YES" if r2v > 0.05 else ("WEAK" if r2v > 0 else "NO")
        print(f"{lbl:<25} {r2v:>8.4f} {sig:>8}")

    # kappa block structure (is it occ-virt dominant for RHF?)
    print(f"\n=== Kappa block structure (RHF, averaged over {n} samples) ===")
    oo, ov, vv = [], [], []
    for s in group:
        kr, ki = s["kr"], s["ki"]
        oo.append(np.linalg.norm(kr[:nocc, :nocc]) + np.linalg.norm(ki[:nocc, :nocc]))
        ov.append(np.linalg.norm(kr[:nocc, nocc:]) + np.linalg.norm(ki[:nocc, nocc:]))
        vv.append(np.linalg.norm(kr[nocc:, nocc:]) + np.linalg.norm(ki[nocc:, nocc:]))
    oo, ov, vv = np.mean(oo), np.mean(ov), np.mean(vv)
    tot = oo + ov + vv
    print(f"  occ-occ:   {oo:.3f} ({100*oo/tot:.1f}%)")
    print(f"  occ-virt:  {ov:.3f} ({100*ov/tot:.1f}%)")
    print(f"  virt-virt: {vv:.3f} ({100*vv/tot:.1f}%)")


if __name__ == "__main__":
    main()
