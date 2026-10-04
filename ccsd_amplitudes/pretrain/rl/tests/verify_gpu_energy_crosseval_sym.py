#!/usr/bin/env python3
"""[verifier] Is the ~5e-10 Ha engine-vs-ffsim difference at norb 16 caused by the float64 asymmetry of the stored
integrals (two_body is 8-fold symmetric only to ~1e-15, one_body symmetric to ~3e-15; every evaluator reads a
different triangle)?  Evaluate one fixed state (engine complex128 LUCJ state, normalized) with the engine, ffsim's
linear operator and pyscf RDMs, for the stored integrals AND for exactly symmetrized integrals."""
from __future__ import annotations
import json, pickle, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402
from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402

name, cand = "C2N2_rxn3923_R", "label"
tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks/baselines.pkl", "rb"))["tasks"]
(_, _, U, Z, t1), = [t for t in tasks if t[0] == (name, cand)]
d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
h1, eri, const = d["one_body"], d["two_body"], float(d["constant"])
norb, k = int(d["norb"]), int(d["nelec_a"])
perms = [(0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 0, 1), (2, 3, 1, 0), (3, 2, 1, 0)]
eri_s = sum(eri.transpose(p) for p in perms) / 8
h1_s = 0.5 * (h1 + h1.T)
print(f"asym: h1 {np.abs(h1 - h1_s).max():.1e}  eri {np.abs(eri - eri_s).max():.1e}", flush=True)
out = {}
eng = LUCJEnergyGPU(h1, eri, const, norb, (k, k), dtype=torch.complex128, max_mem_gb=7.0)
eng.energy(U, Z, t1=t1)
psi = eng.state(U, Z, t1=t1).cpu().numpy()
psi = psi / np.sqrt(np.vdot(psi, psi).real)
out["engine_full"] = eng.energy_of_state(psi)
eng.release(); del eng; torch.cuda.empty_cache()
eng = LUCJEnergyGPU(h1_s, eri_s, const, norb, (k, k), dtype=torch.complex128, max_mem_gb=7.0)
out["engine_sym"] = eng.energy_of_state(psi)
eng.release(); del eng; torch.cuda.empty_cache()
print(json.dumps(out), flush=True)
from pyscf.fci import direct_spin1
t0 = time.time()
dm1 = np.zeros((norb, norb)); dm2 = np.zeros((norb,) * 4)
for part in (psi.real.copy(), psi.imag.copy()):
    a, b = direct_spin1.make_rdm12(part, norb, (k, k))
    dm1 += a; dm2 += b
out["rdm_full"] = const + float(np.einsum("pq,pq->", h1, dm1) + 0.5 * np.einsum("pqrs,pqrs->", eri, dm2))
out["rdm_sym"] = const + float(np.einsum("pq,pq->", h1_s, dm1) + 0.5 * np.einsum("pqrs,pqrs->", eri_s, dm2))
print(f"rdm done ({time.time() - t0:.0f}s)", json.dumps(out), flush=True)
import ffsim
v = psi.reshape(-1)
for tag, (hh, ee) in (("linop_full", (h1, eri)), ("linop_sym", (h1_s, eri_s))):
    t0 = time.time()
    lin = ffsim.linear_operator(ffsim.MolecularHamiltonian(hh, ee, const), norb=norb, nelec=(k, k))
    out[tag] = float(np.vdot(v, lin @ v).real)
    print(f"{tag} ({time.time() - t0:.0f}s): {out[tag]:.12f}", flush=True)
for a in out:
    print(f"  {a:12s} {out[a]:.12f}   - engine_sym {out[a] - out['engine_sym']:+.3e}", flush=True)
json.dump(out, open(ROOT / "pretrain/rl/tests/results/verify_gpu_energy_crosseval_sym.json", "w"), indent=1)
print("done", flush=True)
