#!/usr/bin/env python3
"""[verifier] norb 17 / 18 checks of pretrain/rl/gpu_energy.py on a 12 GB TITAN Xp.

  n17_10_7 : C2HNO (10,7) rl1 tasks: complex64 AND complex128 engine vs the stored full-precision CPU energies
  n17_9_8  : C3HN (9,8) n17_pre tasks (complex64) vs the 1-decimal logged corr%
  n18      : norb 18 (11,7): one physical state through different internal paths (default fused t=6,
             fuse=False, fuse_t=7, a diagonal-phase gauge U_k -> U_k D_k that leaves the LUCJ operator
             unchanged, a different energy tile) + analytic Slater-determinant energies (Z = 0) for new seeds
  oom      : norb 18 (9,9) / (10,8) on 12 GB must raise MemoryError cleanly
Peak memory: torch max_memory_reserved (+ context, see the nvidia-smi sampler log of the driver script).
"""
from __future__ import annotations

import argparse
import json
import pickle
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402

from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac  # noqa: E402


def haar(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def near_id(n, rng, s):
    from scipy.linalg import expm
    a = s * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    return expm(a - a.conj().T)


def slater_energy(h, eri, const, F, k):
    C = F[:, :k]
    D = np.conj(C) @ C.T
    J = np.einsum("pqrs,pq,rs->", eri, D, D)
    K = np.einsum("pqrs,ps,rq->", eri, D, D)
    return float((const + 2 * np.einsum("pq,pq->", h, D) + 2 * J - K).real)


def mem():
    return dict(max_alloc_gb=torch.cuda.max_memory_allocated() / 1e9,
                max_reserved_gb=torch.cuda.max_memory_reserved() / 1e9)


def run_n17_10_7(ntask, out):
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks/n17_rl1.pkl", "rb"))["tasks"]
    ref = json.load(open(ROOT / "pretrain/opt_true/results/energy_n17_rl1.json"))["per_molecule"]
    rng = np.random.default_rng(5)
    pick = [tasks[i] for i in sorted(rng.choice(len(tasks), size=ntask, replace=False))]
    rows = []
    for dt in (torch.complex64, torch.complex128):
        for (name, cand), _, U, Z, t1 in pick:
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
            eng = LUCJEnergyGPU.from_npz(name, dtype=dt)
            t0 = time.time()
            E = eng.energy(U, Z, t1=t1)
            t = time.time() - t0
            r = ref[name][cand]
            row = dict(part="n17_10_7", name=name, cand=cand, dtype=str(dt), E=E, E_ref=float(r["E"]),
                       dE=E - float(r["E"]), dcf=eng.corr_frac(E) - float(r["corr_frac"]), cf=eng.corr_frac(E),
                       t=t, tile=eng.T, mem_plan=eng.mem_plan, **mem())
            rows.append(row)
            print(f"[n17 (10,7)] {name} {cand} {str(dt)[6:]}: E {E:.10f} dE {row['dE']:+.2e} dcf {row['dcf']:+.1e} "
                  f"corr% {100 * row['cf']:.4f}  {t:.1f}s  reserved {row['max_reserved_gb']:.2f} GB", flush=True)
            eng.release()
            del eng
    by = {}
    for r in rows:
        by.setdefault((r["name"], r["cand"]), {})[r["dtype"]] = r["E"]
    for key, d in by.items():
        if len(d) == 2:
            print(f"   c8 - c16 {key}: {d['torch.complex64'] - d['torch.complex128']:+.2e}", flush=True)
    out.extend(rows)


def run_n17_9_8(out):
    logref = {}
    for f in ["rl_runs/energy17_54608122.out", "rl_runs/energy17_54613781.out"]:
        for ln in open(ROOT / f):
            m = re.match(r"\s+(\S+)\s+(\S+)\s+corr%\s+([-\d.]+)", ln)
            if m:
                logref[(m.group(1), m.group(2))] = float(m.group(3)) / 100
    tasks = pickle.load(open(ROOT / "runs_ot/energy_tasks/n17_pre.pkl", "rb"))["tasks"]
    sel = [t for t in tasks if t[1].startswith("C3HN")]
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    eng = None
    for (name, cand), _, U, Z, t1 in sel:
        if eng is None:
            eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex64)
        elif eng.name != name:
            eng.set_npz(name)
        t0 = time.time()
        E = eng.energy(U, Z, t1=t1)
        t = time.time() - t0
        cf = eng.corr_frac(E)
        lr = logref.get((name, cand))
        row = dict(part="n17_9_8", name=name, cand=cand, E=E, cf=cf, cf_log=lr,
                   dcf_log=None if lr is None else cf - lr, t=t, tile=eng.T, k=eng.k, norb=eng.norb,
                   mem_plan=eng.mem_plan, timing=eng.timing, **mem())
        out.append(row)
        print(f"[n17 (9,8)] {name} {cand}: corr% {100 * cf:.4f} (log {'-' if lr is None else f'{100 * lr:.1f}'})"
              f"  {t:.1f}s  reserved {row['max_reserved_gb']:.2f} GB  alloc {row['max_alloc_gb']:.2f} GB", flush=True)
    eng.release()


def run_n18(name, out, seed):
    rng = np.random.default_rng(seed)
    d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
    h, eri, const = d["one_body"], d["two_body"], float(d["constant"])
    n, k = int(d["norb"]), int(d["nelec_a"])
    U = np.stack([near_id(n, rng, 0.25) for _ in range(2)])
    Z = 0.3 * rng.normal(size=(2, n, n))
    Z = 0.5 * (Z + Z.transpose(0, 2, 1))
    t1 = 0.05 * rng.normal(size=(k, n - k))
    Dg = np.exp(1j * rng.uniform(0, 2 * np.pi, size=(2, n)))
    Ug = U * Dg[:, None, :]                        # U_k D_k: same UCJ operator
    variants = [("default", dict(), U), ("gauge", dict(), Ug), ("fuse_t7", dict(fuse_t=7), U),
                ("tile384", dict(tile=384), U), ("nofuse", dict(fuse=False), U)]
    res = {}
    for tag, kw, UU in variants:
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex64, **kw)
        t0 = time.time()
        E = eng.energy(UU, Z, t1=t1)
        t = time.time() - t0
        res[tag] = dict(E=E, t=t, t_top=eng.t_top, tile=eng.T, norm=eng.last["norm"], timing=eng.timing, **mem())
        print(f"[n18 {name} ({k},{n - k})] {tag:8s}: E {E:.10f}  corr% {100 * eng.corr_frac(E):.4f}  t_top {eng.t_top}"
              f"  tile {eng.T}  {t:.1f}s  reserved {res[tag]['max_reserved_gb']:.2f} GB", flush=True)
        if tag == "default":
            # analytic Slater checks with the same engine (Z = 0, new seeds), incl. HF
            sl = []
            for i in range(3):
                Us = np.stack([haar(n, rng) for _ in range(2)]) if i else np.stack([np.eye(n)] * 2).astype(complex)
                tt = 0.3 * rng.normal(size=(k, n - k)) if i else None
                Es = eng.energy(Us, np.zeros((2, n, n)), t1=tt)
                from ffsim.variational.util import orbital_rotation_from_t1_amplitudes
                F = orbital_rotation_from_t1_amplitudes(tt) if tt is not None else np.eye(n)
                Ea = slater_energy(h, eri, const, F, k)
                sl.append(dict(case=i, E=Es, E_analytic=Ea, dE=Es - Ea))
                print(f"     slater {i}: E {Es:.10f} analytic {Ea:.10f} dE {Es - Ea:+.2e}"
                      + (f" (npz e_hf {float(d['e_hf']):.10f})" if i == 0 else ""), flush=True)
            res["slater"] = sl
        eng.release()
        del eng
    base = res["default"]["E"]
    for tag in res:
        if tag != "slater":
            res[tag]["dE_vs_default"] = res[tag]["E"] - base
    print("   spread vs default: " + "  ".join(f"{t}: {res[t]['dE_vs_default']:+.1e}" for t in res if t != "slater"),
          flush=True)
    out.append(dict(part="n18", name=name, k=k, norb=n, res=res))


def run_oom(out):
    for name in ["C4H2_rxn9388_P", "C2H6O_rxn2392_P"]:
        try:
            eng = LUCJEnergyGPU.from_npz(name, dtype=torch.complex64)
            out.append(dict(part="oom", name=name, result="constructed (unexpected)", mem_plan=eng.mem_plan))
            eng.release()
        except MemoryError as e:
            out.append(dict(part="oom", name=name, result=f"MemoryError: {e}"))
            print(f"[oom] {name}: MemoryError: {e}", flush=True)
        torch.cuda.empty_cache()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["n17_10_7", "n17_9_8", "n18", "oom"])
    ap.add_argument("--n17-tasks", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = []
    if "oom" in args.parts:
        run_oom(out)
    if "n17_10_7" in args.parts:
        run_n17_10_7(args.n17_tasks, out)
    if "n17_9_8" in args.parts:
        run_n17_9_8(out)
    if "n18" in args.parts:
        run_n18("C2H2O2_rxn3279_TS", out, seed=31)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=1, default=str)
    print("->", args.out)


if __name__ == "__main__":
    main()
