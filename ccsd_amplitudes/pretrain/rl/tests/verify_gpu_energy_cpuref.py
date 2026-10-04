#!/usr/bin/env python3
"""[verifier] Fresh CPU ffsim reference for stored energy tasks (same code path as the CPU energy job:
SVD polar of U, exact_energy(make_ucj_op(Z, U, "square", t1))).  CPU only; pin threads via env.

  OMP_NUM_THREADS=4 ... python -m pretrain.rl.tests.verify_gpu_energy_cpuref --pkl baselines \
      --name C2N2_rxn3923_R --cand label --out <json>
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pretrain.rl.energy import make_ucj_op  # noqa: E402
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--cand", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import ffsim
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{args.pkl}.pkl", "rb"))["tasks"]
    sel = [t for t in tasks if t[0] == (args.name, args.cand)]
    assert len(sel) == 1, len(sel)
    _, name, U, Z, t1 = sel[0]
    ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
    W, _, Vh = np.linalg.svd(U.astype(np.complex128))
    U = W @ Vh
    t0 = time.time()
    op = make_ucj_op(Z, U, "square", t1=t1)
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, nelec), op, norb=norb, nelec=nelec)
    t_state = time.time() - t0
    linop = ffsim.linear_operator(ham, norb=norb, nelec=nelec)
    hpsi = linop @ psi
    E = float(np.vdot(psi, hpsi).real)
    nrm = float(np.vdot(psi, psi).real)
    t_all = time.time() - t0
    out = dict(pkl=args.pkl, name=name, cand=args.cand, norb=norb, nelec=list(nelec), E=E, norm=nrm,
               E_normalized=E / nrm, corr_frac=(e_hf - E) / (e_hf - e_ccsd), t_state=t_state, t=t_all)
    print(json.dumps(out), flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
