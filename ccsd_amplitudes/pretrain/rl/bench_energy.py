#!/usr/bin/env python3
"""Benchmark the two energy evaluators on cached Hamiltonians.

For each molecule: builds the exact-DF (optimize=False) and compressed
(optimize=True) LUCJ operators for the given connectivity, then reports
  - exact statevector energy (only if norb <= --exact-max-norb)
  - MPS+SQD energy for a grid of (max_bond, shots)
with wall times.  Energies are given as % of CCSD correlation energy recovered.

Usage (from ccsd_amplitudes/, ffsim venv):
    python -m pretrain.rl.bench_energy --names C2H3N_rxn2857_TS C2H5N3O_rxn3552_P
"""
from __future__ import annotations

import argparse
import os
import time

os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np  # noqa: E402

from pretrain.rl.energy import exact_energy, make_ucj_op, mps_sqd_energy, z_from_op  # noqa: E402
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--ham-dir", default="rhf_hamiltonians")
    ap.add_argument("--connectivity", default="square")
    ap.add_argument("--exact-max-norb", type=int, default=16)
    ap.add_argument("--max-bonds", nargs="+", type=int, default=[32, 64, 128])
    ap.add_argument("--shots", nargs="+", type=int, default=[1000, 3000])
    ap.add_argument("--no-t1", action="store_true")
    args = ap.parse_args()

    import ffsim

    for name in args.names:
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(args.ham_dir, name)
        d = np.load(f"{args.data_dir}/{name}.npz")
        t2 = d["t2"].astype(float)
        t1 = None if args.no_t1 else d["t1"].astype(float)
        ecorr = e_ccsd - e_hf
        print(f"\n=== {name}: norb={norb} nelec={nelec} E_HF={e_hf:.6f} "
              f"E_CCSD={e_ccsd:.6f} (corr {ecorr:.5f})")

        from pretrain.rl.energy import interaction_pairs
        pairs = interaction_pairs(args.connectivity, norb)
        ops = {}
        t0 = time.time()
        op0 = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
            t2, t1=t1, n_reps=2, interaction_pairs=pairs, optimize=False)
        ops["exact-DF"] = op0
        print(f"  exact DF built in {time.time()-t0:.1f}s")
        t0 = time.time()
        op1 = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
            t2, t1=t1, n_reps=2, interaction_pairs=pairs, optimize=True,
            options=dict(maxiter=300))
        ops["compressed"] = op1
        print(f"  compressed DF built in {time.time()-t0:.1f}s")
        # round trip through (Z, U) -> make_ucj_op to validate that helper.
        # ffsim splits the optimized (tridiagonal) Z into the alpha-alpha block
        # (off-diagonal pairs) and the alpha-beta block (diagonal pairs); the
        # label pipeline stores the full Z, so recombine the disjoint blocks.
        Z_full = z_from_op(op0, args.connectivity)
        op0b = make_ucj_op(Z_full, op0.orbital_rotations, args.connectivity, t1=t1)
        assert np.allclose(op0b.diag_coulomb_mats, op0.diag_coulomb_mats), "make_ucj_op mismatch"

        for label, op in ops.items():
            if norb <= args.exact_max_norb:
                t0 = time.time()
                e = exact_energy(ham, norb, nelec, op)
                print(f"  [{label:10s}] exact      E={e:.6f}  corr%={(e_hf-e)/(-ecorr)*100:6.1f}"
                      f"  ({time.time()-t0:.1f}s)")
            for mb in args.max_bonds:
                for sh in args.shots:
                    t0 = time.time()
                    try:
                        e, info = mps_sqd_energy(ham, norb, nelec, op, shots=sh,
                                                 max_bond=mb, seed=0)
                        print(f"  [{label:10s}] MPS chi={mb:4d} shots={sh:5d}  E={e:.6f}"
                              f"  corr%={(e_hf-e)/(-ecorr)*100:6.1f}  "
                              f"t_mps={info['t_mps']:.1f} t_smp={info['t_sample']:.1f} "
                              f"t_sqd={info['t_sqd']:.1f} uniq={info['n_unique']} "
                              f"chi_max={info['max_bond_reached']}  total {time.time()-t0:.1f}s",
                              flush=True)
                    except Exception as ex:  # noqa: BLE001
                        print(f"  [{label}] MPS chi={mb} shots={sh} FAILED: "
                              f"{type(ex).__name__}: {ex}", flush=True)


if __name__ == "__main__":
    main()
