#!/usr/bin/env python3
"""Task 5: accuracy and cost of the density-matrix TN engine (pretrain/rl/tn_energy.py, method "dm") beyond norb 29.

Runs a list of (molecule, candidate, chi) jobs in order in ONE process (one engine per molecule, reused across chi, so
the block2 MPO is built once per molecule and its build time is recorded on the first job of that molecule).  One
json line per job with E, % of the CCSD correlation energy, truncation weights, every timing of the engine
(t_state = MPS build incl. t_env / t_sweep, t_mpo_build, t_convert, t_expect), torch peak allocated / reserved memory,
wall-clock start / end (to match the nvidia-smi footprint log of pretrain/followups/gpu_monitor.py), host load, and
every engine setting.

Parameters: a policy_dump pickle (--params, tasks ((name, tag), name, U, Z, t1)) and/or n29 bench parameter files
(--bench-params name:tag:file.npz with U_<tag>, Z_<tag>; t1 from rhf_dataset).

Usage (thread caps must be set in the environment before the start, e.g. OMP_NUM_THREADS=4 ...):
  CUDA_VISIBLE_DEVICES=4 python3 pretrain/followups/task5_dm_bench.py --params P.pkl \
      --jobs C7H16_rxn8743_R:rl4n29f:64,C7H16_rxn8743_R:rl4n29f:128 --b2-threads 4 --mem-cap-gb 9 --out out.jsonl
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    assert _v in os.environ, f"set {_v} explicitly (CPU budget)"
import argparse  # noqa: E402
import gc  # noqa: E402
import importlib  # noqa: E402
import json  # noqa: E402
import pickle  # noqa: E402
import socket  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl.hamiltonian import load_hamiltonian  # noqa: E402


def load_params(args):
    P = {}
    if args.params:
        d = pickle.load(open(ROOT / args.params, "rb"))
        for (key, name, U, Z, t1) in d["tasks"]:
            P[(name, key[1])] = (U, Z, t1)
    for spec in args.bench_params or []:
        name, tag, f = spec.split(":")
        B = np.load(ROOT / f)
        t1 = np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t1"].astype(np.float64)
        P[(name, tag)] = (B[f"U_{tag}"], B[f"Z_{tag}"], t1)
    return P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=None)
    ap.add_argument("--bench-params", action="append", default=None)
    ap.add_argument("--jobs", required=True, help="name:tag:chi,... (run in this order)")
    ap.add_argument("--impl", default="current", choices=["current", "v1"])
    ap.add_argument("--method", default=None, help="dm / zipup (default: engine default = dm on CUDA)")
    ap.add_argument("--b2-threads", type=int, default=4)
    ap.add_argument("--stack-mem-gb", type=float, default=2.0)
    ap.add_argument("--mem-cap-gb", type=float, default=None, help="torch per-process memory cap (reserved)")
    ap.add_argument("--basis-cache", default=str(ROOT / "rl_runs/followups/task5_dm_scaling/basis_cache"))
    ap.add_argument("--label", default="", help="free-form run label stored with every record")
    ap.add_argument("--monitor-log", default=None, help="gpu_monitor.py log of this GPU: repeat a job (up to "
                    "--max-tries) while the whole GPU was >= 95 %% busy for more than --busy-max of its state build "
                    "(another user's job shares the GPU)")
    ap.add_argument("--busy-max", type=float, default=0.2)
    ap.add_argument("--max-tries", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    T = importlib.import_module("pretrain.rl.tn_energy_v1" if args.impl == "v1" else "pretrain.rl.tn_energy")
    P = load_params(args)
    jobs = []
    for spec in args.jobs.split(","):
        name, tag, chi = spec.split(":")
        assert (name, tag) in P, f"no parameters for {name}:{tag}"
        jobs.append((name, tag, int(chi)))
    dev_total = torch.cuda.get_device_properties(0).total_memory
    if args.mem_cap_gb:
        torch.cuda.set_per_process_memory_fraction(min(1.0, args.mem_cap_gb * (1 << 30) / dev_total), 0)
    gpu = torch.cuda.get_device_name(0)
    scratch_root = Path(tempfile.gettempdir()) / f"task5_{os.getpid()}"
    scratch_root.mkdir(parents=True, exist_ok=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cur, ev, t_init = None, None, 0.0

    def busy_frac(t0, t1):
        if not args.monitor_log:
            return None
        time.sleep(2.5)                                    # let the monitor write the last samples
        u = []
        for ln in open(args.monitor_log):
            try:
                m = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if t0 <= m.get("t", 0) <= t1 and "gpu_util" in m:
                u.append(m["gpu_util"])
        return float(np.mean(np.array(u) >= 95)) if u else None

    jobs = [j + (k,) for j in jobs for k in range(1)]
    queue = list(jobs)
    tries = {}
    while queue:
        name, tag, chi, _ = queue.pop(0)
        tries[(name, tag, chi)] = tries.get((name, tag, chi), 0) + 1
        ham, norb, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
        rec = {"name": name, "cand": tag, "norb": norb, "nelec": list(nelec), "chi": chi, "e_hf": e_hf,
               "e_ccsd": e_ccsd, "pid": os.getpid(), "host": socket.gethostname(), "gpu": gpu,
               "cuda_visible": os.environ.get("CUDA_VISIBLE_DEVICES"), "omp": os.environ["OMP_NUM_THREADS"],
               "impl": args.impl, "mem_cap_gb": args.mem_cap_gb, "label": args.label}
        try:
            if cur != name:
                ev = None
                gc.collect()
                torch.cuda.empty_cache()
                kw = {} if args.method is None else {"method": args.method}
                t0 = time.time()
                ev = T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec,
                                    max_bond=chi, device="cuda", name=name, block2_threads=args.b2_threads,
                                    basis_cache=args.basis_cache, stack_mem=int(args.stack_mem_gb * (1 << 30)),
                                    scratch=tempfile.mkdtemp(prefix="b2_", dir=scratch_root), **kw)
                t_init = time.time() - t0
                cur = name
            U, Z, t1 = P[(name, tag)]
            torch.cuda.reset_peak_memory_stats()
            load0 = os.getloadavg()[0]
            t_start = time.time()
            E, info = ev.energy(U, Z, t1, max_bond=chi)
            t_end = time.time()
            bf = busy_frac(t_start, t_start + info["t_state"])
            rec["busy_frac"], rec["attempt"] = bf, tries[(name, tag, chi)]
            if bf is not None and bf > args.busy_max and tries[(name, tag, chi)] < args.max_tries:
                queue.insert(0, (name, tag, chi, 0))         # repeat right away (engine and MPO are cached)
            rec.update({"E": E, "corr_pct": (e_hf - E) / (e_hf - e_ccsd) * 100, "t_init": t_init,
                        "t_start": t_start, "t_end": t_end, "wall": t_end - t_start, "load_start": load0,
                        "load_end": os.getloadavg()[0],
                        "torch_peak_alloc_MB": torch.cuda.max_memory_allocated() / 2 ** 20,
                        "torch_peak_reserved_MB": torch.cuda.max_memory_reserved() / 2 ** 20,
                        **ev.settings(),
                        **{k: v for k, v in info.items() if not isinstance(v, (list, dict))}})
            t_init = 0.0
        except Exception as e:  # noqa: BLE001
            rec.update({"error": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()[-2000:],
                        "t_end": time.time()})
            cur, ev = None, None                      # rebuild the engine after a failure
            gc.collect()
            torch.cuda.empty_cache()
        print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in rec.items()
                          if k in ("name", "cand", "chi", "E", "corr_pct", "discarded_sum", "t_init", "t_state",
                                   "t_mpo_build", "t_expect", "wall", "torch_peak_alloc_MB",
                                   "torch_peak_reserved_MB", "busy_frac", "attempt", "error")}), flush=True)
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
