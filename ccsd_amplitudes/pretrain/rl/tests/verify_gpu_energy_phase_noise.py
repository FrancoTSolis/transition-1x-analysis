#!/usr/bin/env python3
"""[verifier] Float64 noise of the evaluators at norb 16: E(e^{i theta} psi) = E(psi) exactly (H real symmetric), but the
real/imaginary split (which pyscf contracts separately) and every rounding pattern changes with theta.  The spread over
theta is a direct measure of each evaluator's rounding noise (engine complex128 vs ffsim linear operator)."""
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
eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex128, max_mem_gb=7.0)
eng.energy(U, Z, t1=t1)
psi = eng.state(U, Z, t1=t1).cpu().numpy()
psi = psi / np.sqrt(np.vdot(psi, psi).real)
thetas = [0.0, np.pi / 4, np.pi / 3, 1.0, 2.5]
out = dict(engine={}, linop={})
for th in thetas:
    out["engine"][f"{th:.4f}"] = eng.energy_of_state(psi * np.exp(1j * th))
print("engine:", json.dumps(out["engine"]), " spread %.2e" % (max(out["engine"].values()) - min(out["engine"].values())), flush=True)
eng.release(); del eng; torch.cuda.empty_cache()
import ffsim
d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
ham = ffsim.MolecularHamiltonian(d["one_body"], d["two_body"], float(d["constant"]))
lin = ffsim.linear_operator(ham, norb=int(d["norb"]), nelec=(int(d["nelec_a"]),) * 2)
for th in (np.pi / 4, 1.0):
    t0 = time.time()
    v = (psi * np.exp(1j * th)).reshape(-1)
    out["linop"][f"{th:.4f}"] = float(np.vdot(v, lin @ v).real)
    print(f"linop theta={th:.4f}: {out['linop'][f'{th:.4f}']:.12f}  ({time.time() - t0:.0f}s)", flush=True)
out["linop"]["0.0000 (earlier run)"] = -182.394146072784
json.dump(out, open(ROOT / "pretrain/rl/tests/results/verify_gpu_energy_phase_noise.json", "w"), indent=1)
print("done", flush=True)
