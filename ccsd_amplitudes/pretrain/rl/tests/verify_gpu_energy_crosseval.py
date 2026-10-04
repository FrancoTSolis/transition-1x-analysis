#!/usr/bin/env python3
"""[verifier] Cross-evaluate one state with three float64 evaluators to arbitrate small (1e-10) differences:
engine c16 evaluator, ffsim linear operator <psi|H|psi>, and pyscf 1-/2-RDM contraction (real and imaginary
parts of the complex CI vector separately; H is real symmetric so the cross terms cancel).
The state is the engine's complex128 LUCJ state of a stored task (normalized)."""
from __future__ import annotations
import argparse, json, pickle, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402
from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", default="baselines"); ap.add_argument("--name", default="C2N2_rxn3923_R")
    ap.add_argument("--cand", default="label"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks" / f"{a.pkl}.pkl", "rb"))["tasks"]
    (_, _, U, Z, t1), = [t for t in tasks if t[0] == (a.name, a.cand)]
    eng = LUCJEnergyGPU.from_npz(a.name, dtype=torch.complex128, max_mem_gb=7.0)
    E_eng = eng.energy(U, Z, t1=t1)
    psi = eng.state(U, Z, t1=t1).cpu().numpy()
    asym = float(np.abs(psi - psi.T).max())
    nrm = float(np.vdot(psi, psi).real)
    psi = psi / np.sqrt(nrm)
    E_eng_state = eng.energy_of_state(psi)
    eng.release(); del eng; torch.cuda.empty_cache()
    print(f"engine: E {E_eng:.12f}  energy_of_state(normalized) {E_eng_state:.12f}  |psi|^2-1 {nrm-1:+.1e}  asym {asym:.1e}", flush=True)
    import ffsim
    from pyscf.fci import direct_spin1
    d = np.load(ROOT / "rhf_hamiltonians" / f"{a.name}.npz")
    h1, eri, const = d["one_body"], d["two_body"], float(d["constant"])
    norb, k = int(d["norb"]), int(d["nelec_a"])
    ham = ffsim.MolecularHamiltonian(h1, eri, const)
    t0 = time.time()
    lin = ffsim.linear_operator(ham, norb=norb, nelec=(k, k))
    v = psi.reshape(-1)
    E_lin = float(np.vdot(v, lin @ v).real)
    t_lin = time.time() - t0
    print(f"ffsim linop: E {E_lin:.12f}  ({t_lin:.0f}s)  engine - linop {E_eng_state - E_lin:+.3e}", flush=True)
    t0 = time.time()
    E_rdm = const
    for part in (psi.real.copy(), psi.imag.copy()):
        dm1, dm2 = direct_spin1.make_rdm12(part, norb, (k, k))
        E_rdm += float(np.einsum("pq,pq->", h1, dm1) + 0.5 * np.einsum("pqrs,pqrs->", eri, dm2))
    t_rdm = time.time() - t0
    print(f"pyscf rdm12: E {E_rdm:.12f}  ({t_rdm:.0f}s)  engine - rdm {E_eng_state - E_rdm:+.3e}  linop - rdm {E_lin - E_rdm:+.3e}", flush=True)
    out = dict(name=a.name, cand=a.cand, E_engine=E_eng, E_engine_state=E_eng_state, E_linop=E_lin, E_rdm=E_rdm,
               norm_minus_1=nrm - 1, asym=asym, t_lin=t_lin, t_rdm=t_rdm)
    json.dump(out, open(a.out, "w"), indent=1)
    print("->", a.out)

if __name__ == "__main__":
    main()
