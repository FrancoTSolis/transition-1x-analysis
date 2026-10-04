#!/usr/bin/env python3
"""Random-case validation of pretrain/rl/gpu_energy.py against ffsim (exact_energy of make_ucj_op).

Cases: norb 4..10, nalpha = nbeta = k with k != norb - k included, random Haar U (2 reps), random real
symmetric Z, random t1, Hamiltonians = random real tensors with 8-fold symmetry (indefinite ERI
matrix: exercises negative eigenvalues) or random PSD ERIs, plus real pyscf molecules (STO-3G).
Connectivity: mostly "square", some "all-to-all" / "hex".

  CUDA_VISIBLE_DEVICES=2 OMP_NUM_THREADS=4 python -m pretrain.rl.tests.test_gpu_energy_random \
      --n-cases 36 --out pretrain/rl/tests/results/gpu_energy_random.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import torch  # noqa: E402

from pretrain.rl.energy import exact_energy, make_ucj_op  # noqa: E402
from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402


def sym8(t):
    perms = [(0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2),
             (2, 3, 0, 1), (3, 2, 0, 1), (2, 3, 1, 0), (3, 2, 1, 0)]
    return sum(t.transpose(p) for p in perms) / 8


def random_ham(norb, rng, kind):
    h = rng.normal(size=(norb, norb))
    h = 0.5 * (h + h.T)
    if kind == "sym8":
        eri = sym8(rng.normal(size=(norb,) * 4)) * 0.5
    else:  # PSD (Coulomb-like)
        L = rng.normal(size=(norb * norb, norb, norb))
        L = 0.5 * (L + L.transpose(0, 2, 1))
        eri = np.einsum("lpq,lrs->pqrs", L, L) / (norb * norb)
    return h, eri, float(rng.normal())


def pyscf_ham(name, norb_max):
    import ffsim
    from pyscf import gto, scf
    geoms = {
        "H6": "".join(f"H 0 0 {1.0 * i}\n" for i in range(6)),
        "H8": "".join(f"H 0 0 {0.9 * i}\n" for i in range(8)),
        "H10": "".join(f"H 0 0 {0.95 * i}\n" for i in range(10)),
        "LiH": "Li 0 0 0\nH 0 0 1.6",
        "BeH2": "Be 0 0 0\nH 0 0 1.3\nH 0 0 -1.3",
        "H2O": "O 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59",
    }
    mol = gto.M(atom=geoms[name], basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    data = ffsim.MolecularData.from_scf(mf)
    ham = data.hamiltonian
    return (ham.one_body_tensor.real, ham.two_body_tensor.real, ham.constant, data.norb, data.nelec)


def haar(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cases", type=int, default=36)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--backends", nargs="+", default=["numba", "torch"])
    ap.add_argument("--out", default="pretrain/rl/tests/results/gpu_energy_random.json")
    args = ap.parse_args()
    import ffsim

    rng = np.random.default_rng(args.seed)
    shapes = [(4, 1), (4, 2), (5, 2), (5, 3), (6, 2), (6, 3), (6, 4), (7, 2), (7, 3), (7, 4), (7, 5),
              (8, 3), (8, 4), (8, 5), (8, 6), (9, 3), (9, 4), (9, 5), (9, 6), (10, 3), (10, 4), (10, 5),
              (10, 6), (10, 7)]
    cases = []
    for c in range(args.n_cases):
        norb, k = shapes[c % len(shapes)]
        kind = "sym8" if c % 2 == 0 else "psd"
        conn = ["square", "square", "square", "all-to-all", "hex"][c % 5]
        cases.append(("random", norb, k, kind, conn))
    for mol in ["H6", "H8", "LiH", "BeH2", "H2O", "H10"]:
        cases.append(("pyscf", mol, None, None, "square"))

    rows = []
    worst = {}
    for ci, case in enumerate(cases):
        if case[0] == "random":
            _, norb, k, kind, conn = case
            h, eri, const = random_ham(norb, rng, kind)
            label = f"rand n={norb} k={k} {kind} {conn}"
        else:
            h, eri, const, norb, nelec = pyscf_ham(case[1], 10)
            k = nelec[0]
            conn = case[4]
            label = f"pyscf {case[1]} n={norb} k={k} {conn}"
        U = np.stack([haar(norb, rng) for _ in range(2)])
        Z = rng.normal(size=(2, norb, norb))
        Z = 0.5 * (Z + Z.transpose(0, 2, 1))
        t1 = 0.2 * rng.normal(size=(k, norb - k)) if norb - k > 0 else None
        ham = ffsim.MolecularHamiltonian(h, eri, const)
        op = make_ucj_op(Z, U, conn, t1=t1)
        t0 = time.time()
        E_ref = exact_energy(ham, norb, (k, k), op)
        psi_ref = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, (k, k)), op, norb=norb, nelec=(k, k))
        t_ref = time.time() - t0
        row = dict(case=label, norb=norb, k=k, conn=conn, E_ref=E_ref, t_ref=t_ref)
        for be in args.backends:
            for dt in (torch.complex128, torch.complex64):
                tag = f"{be}_{'c16' if dt == torch.complex128 else 'c8'}"
                eng = LUCJEnergyGPU(h, eri, const, norb, (k, k), dtype=dt, backend=be, max_mem_gb=2.0)
                E = eng.energy(U, Z, t1=t1, connectivity=conn)
                psi = eng.state(U, Z, t1=t1, connectivity=conn).cpu().numpy().reshape(-1)
                dpsi = float(np.abs(psi - psi_ref).max())
                row[tag] = dict(E=E, dE=E - E_ref, dpsi=dpsi, nq=eng.nq, npos=eng.npos)
                w = worst.setdefault(tag, dict(dE=0.0, dpsi=0.0))
                w["dE"] = max(w["dE"], abs(E - E_ref))
                w["dpsi"] = max(w["dpsi"], dpsi)
                eng.release()
        rows.append(row)
        msg = "  ".join(f"{t}: dE {row[t]['dE']:+.1e} dpsi {row[t]['dpsi']:.1e}" for t in row if "_c" in t)
        print(f"[{ci:2d}] {label:40s} E_ref {E_ref:+.10f}  {msg}", flush=True)
    print("\n=== worst |dE| (Hartree) / max |dpsi| over", len(rows), "cases ===")
    for t, w in worst.items():
        print(f"  {t:12s} max|dE| {w['dE']:.2e}   max|dpsi| {w['dpsi']:.2e}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(worst=worst, rows=rows, args=vars(args)), open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
