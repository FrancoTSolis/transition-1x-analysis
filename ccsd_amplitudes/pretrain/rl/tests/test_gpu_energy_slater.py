#!/usr/bin/env python3
"""Full-size exactness check of pretrain/rl/gpu_energy.py without any CPU statevector (works for norb 17/18).

With Z = 0 the UCJ operator collapses: the engine still applies the dense rotations W_1 = U_1^dag,
W_2 = U_2^dag U_1, W_final = F U_2 one after the other (Slater-determinant init, fused Givens passes,
transposes, energy tiles), but the product is F = exp(t1 - t1^dag), so the state is the closed-shell Slater
determinant with occupied orbitals F[:, :k] and its energy is analytic:
    D_pq = sum_{i<k} conj(F_pi) F_qi,   E = const + 2 sum h_pq D_pq + 2 J - K,
    J = sum (pq|rs) D_pq D_rs,   K = sum (pq|rs) D_ps D_rq.
With t1 = None it must give E_HF of the npz (checks the integrals as well).

  CUDA_VISIBLE_DEVICES=2 python -m pretrain.rl.tests.test_gpu_energy_slater --names C3HN_rxn2586_P \
      --out pretrain/rl/tests/results/gpu_energy_slater.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402

from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402
from pretrain.rl.tests.bench_gpu_energy import haar  # noqa: E402


def slater_energy(h, eri, const, F, k):
    C = F[:, :k]
    D = np.conj(C) @ C.T                       # D_pq = sum_i conj(F_pi) F_qi = <a^dag_p a_q> (one spin)
    J = np.einsum("pqrs,pq,rs->", eri, D, D)
    K = np.einsum("pqrs,ps,rq->", eri, D, D)
    return float((const + 2 * np.einsum("pq,pq->", h, D) + 2 * J - K).real)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--n-random", type=int, default=2)
    ap.add_argument("--t1-scale", type=float, default=0.3)
    ap.add_argument("--dtype", default="c8")
    ap.add_argument("--max-mem-gb", type=float, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
    dt = torch.complex64 if args.dtype == "c8" else torch.complex128
    rng = np.random.default_rng(11)
    rows = []
    for name in args.names:
        d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
        h, eri, const = d["one_body"], d["two_body"], float(d["constant"])
        torch.cuda.reset_peak_memory_stats()
        eng = LUCJEnergyGPU.from_npz(name, dtype=dt, max_mem_gb=args.max_mem_gb)
        n, k = eng.norb, eng.k
        cases = [("hf", None, np.stack([np.eye(n)] * 2).astype(complex))]
        for i in range(args.n_random):
            cases.append((f"rand{i}", args.t1_scale * rng.normal(size=(k, n - k)),
                          np.stack([haar(n, rng) for _ in range(2)])))
        for tag, t1, U in cases:
            Z = np.zeros((2, n, n))
            t0 = time.time()
            E = eng.energy(U, Z, t1=t1)
            dt_s = time.time() - t0
            F = orbital_rotation_from_t1_amplitudes(t1) if t1 is not None else np.eye(n)
            E_ref = slater_energy(h, eri, const, F, k)
            row = dict(name=name, norb=n, k=k, dim=eng.dim, case=tag, E=E, E_ref=E_ref, dE=E - E_ref,
                       e_hf=float(d["e_hf"]), t=dt_s, peak_gb=torch.cuda.max_memory_allocated() / 1e9,
                       dtype=args.dtype, norm=eng.last["norm"])
            rows.append(row)
            print(f"{name} n={n} k={k} dim={eng.dim} {tag:6s}: E {E:.10f}  analytic {E_ref:.10f}  dE {E - E_ref:+.2e}"
                  + (f"  (E_HF npz {d['e_hf']:.10f}, d {E_ref - d['e_hf']:+.1e})" if tag == "hf" else "")
                  + f"  |psi|^2-1 {eng.last['norm'] - 1:+.1e}  {dt_s:.1f}s  peak {row['peak_gb']:.2f} GB", flush=True)
        eng.release()
        del eng
        torch.cuda.empty_cache()
    print("max |dE| =", max(abs(r["dE"]) for r in rows))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(rows, open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
