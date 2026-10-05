#!/usr/bin/env python3
"""Objective server for task 2: one molecule, one GPU engine, requests over stdin/stdout (length-prefixed pickles).

Runs in its own process so that NOMAD (whose PSD-MADS subproblem threads must run with OMP_NUM_THREADS=1: with
more OpenMP threads PyNomad 4.5.1 crashed on a toy problem) never shares OpenMP / CUDA state with torch, numba
and pyscf.  The driver serializes the requests.  Selected-CI diagonalizations run in a further child process
without torch (task2_sci_worker; see the note there on wrong PySCF SCI energies in mixed processes).

Requests: {"cmd": "lucj", "x"}            -> exact LUCJ energy (GPU, complex64)
          {"cmd": "qsci", "x", "shots", "sample_seed", "settings", "verify"}
                                           -> QSCI energy on `shots` exact samples (fixed sampling seed, as Lin
                                              et al.'s quimb sample(seed=0)); the subsampling rng is the SCI
                                              worker's persistent rng.  verify: the full-space Rayleigh quotient
                                              of the returned SCI eigenvector on the GPU engine (complex64),
                                              an independent check of the subspace energy
          {"cmd": "marginals", "x"}         -> exact LUCJ energy and exact single-spin string marginals
          {"cmd": "state_ci", "x", "sampler"}  (scoring) -> exact LUCJ energy, entropy, p_HF and Task-1 CI strings
          {"cmd": "solve", "ci", "spin_sq"}  (scoring) -> SCI energies of the given subspaces (SCI child) + the
                                              full-space check of the lowest one
          {"cmd": "verify_vec", "vec"}       -> full-space <v|H|v> of an SCI eigenvector
          {"cmd": "quit"}
"""
from __future__ import annotations

import argparse
import os
import pickle
import struct
import subprocess
import sys
import time
import traceback
from pathlib import Path


def _read(f):
    h = f.read(8)
    if len(h) < 8:
        return None
    (n,) = struct.unpack("<Q", h)
    return pickle.loads(f.read(n))


def _write(f, obj):
    b = pickle.dumps(obj)
    f.write(struct.pack("<Q", len(b)) + b)
    f.flush()


class SCIChild:
    def __init__(self, name, threads, entropy, log):
        root = Path(__file__).resolve().parents[2]
        env = dict(os.environ)
        for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS",
                  "NUMBA_NUM_THREADS"):
            env[v] = str(threads)
        self.p = subprocess.Popen([sys.executable, "-m", "pretrain.followups.task2_sci_worker", "--name", name,
                                   "--threads", str(threads), "--entropy", str(entropy)],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, cwd=str(root), env=env)
        info = _read(self.p.stdout)
        if not info or not info.get("ready"):
            raise RuntimeError(f"SCI worker failed to start: {info}")

    def call(self, req):
        _write(self.p.stdin, req)
        r = _read(self.p.stdout)
        if r is None:
            raise RuntimeError("SCI worker died")
        if "error" in r:
            raise RuntimeError("SCI worker: " + r["error"] + "\n" + r.get("tb", ""))
        return r

    def close(self):
        try:
            _write(self.p.stdin, {"cmd": "quit"})
            self.p.wait(timeout=60)
        except Exception:  # noqa: BLE001
            self.p.kill()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--max-mem-gb", type=float, default=5.0)
    ap.add_argument("--entropy", type=int, default=0)
    ap.add_argument("--no-sci", action="store_true", help="LUCJ energies only (no SCI child)")
    args = ap.parse_args()
    # protocol on the original stdout; everything printed by libraries goes to stderr
    proto_out = os.fdopen(os.dup(1), "wb")
    os.dup2(2, 1)
    proto_in = sys.stdin.buffer

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import numpy as np
    import torch
    torch.set_num_threads(args.threads)
    from pyscf.fci import cistring
    from pretrain.followups import task2_common as C
    from pretrain.rl.gpu_energy import LUCJEnergyGPU

    ham = C.load_ham(args.name)
    t1 = C.load_t1(args.name)
    eng = LUCJEnergyGPU.from_npz(args.name, max_mem_gb=args.max_mem_gb)
    sci = None if args.no_sci else SCIChild(args.name, args.threads, args.entropy, sys.stderr)
    norb, k = ham["norb"], ham["nelec"][0]

    def full_space_energy(vec):
        """<v|H|v>/<v|v> of an SCI eigenvector embedded in the full CI matrix (engine buffer, complex64)."""
        amp, sa, sb = vec
        ia = torch.as_tensor(cistring.strs2addr(norb, k, np.asarray(sa, dtype=np.int64)), device=eng.device)
        ib = torch.as_tensor(cistring.strs2addr(norb, k, np.asarray(sb, dtype=np.int64)), device=eng.device)
        with torch.no_grad(), torch.cuda.device(eng.device):
            eng.psi.zero_()
            eng.psi[ia[:, None], ib[None, :]] = torch.as_tensor(np.asarray(amp), device=eng.device).to(eng.cdtype)
            asym = float((eng.psi - eng.psi.T).abs().max())
            return float(eng._energy_of_psi()), asym

    _write(proto_out, {"ready": True, "norb": norb, "nelec": ham["nelec"], "e_hf": ham["e_hf"],
                       "e_ccsd": ham["e_ccsd"], "gpu": torch.cuda.get_device_name(), "mem_plan": eng.mem_plan})
    while True:
        req = _read(proto_in)
        if req is None or req.get("cmd") == "quit":
            break
        t0 = time.time()
        try:
            if req["cmd"] == "verify_vec":
                E, asym = full_space_energy(req["vec"])
                out = {"E_full": E, "asym": asym, "t": time.time() - t0}
                _write(proto_out, out)
                continue
            if req["cmd"] == "solve":                # scoring: SCI of given subspaces in the SCI child + check
                r = sci.call({"cmd": "solve", "ci": req["ci"], "spin_sq": req.get("spin_sq"), "return_vec": True})
                vec = r.pop("vec", None)
                if vec is not None:
                    E_full, asym = full_space_energy(vec)
                    r.update(E_full=E_full, d_full=E_full - r["E_vec"], asym=asym)
                r["t"] = time.time() - t0
                _write(proto_out, r)
                continue
            U, Z = C.x_to_uz(req["x"], norb)
            if req["cmd"] == "lucj":
                E = eng.energy(U, Z, t1=t1)
                out = {"f": float(E), "t": time.time() - t0}
            elif req["cmd"] == "qsci":
                psi = eng.state(U, Z, t1=t1)
                bs = C.sample_ci_matrix(psi, norb, k, int(req["shots"]), seed=int(req["sample_seed"]))
                t_samp = time.time() - t0
                q = sci.call({"cmd": "qsci", "bitstrings": bs, "settings": req["settings"],
                              "return_vec": bool(req.get("verify", True))})
                vec = q.pop("vec", None)
                out = {"f": q["E"], **{kk: v for kk, v in q.items() if kk != "E"}, "t_sample": t_samp}
                if vec is not None:
                    E_full, asym = full_space_energy(vec)
                    out.update(E_full=E_full, d_full=E_full - q["E"], asym=asym)
                out["t"] = time.time() - t0
            elif req["cmd"] == "marginals":
                # exact single-spin string marginals p(a) = sum_b |psi[a, b]|^2 (psi = psi^T: alpha = beta)
                E = eng.energy(U, Z, t1=t1)
                with torch.no_grad():
                    pa = torch.zeros(eng.dim, device=eng.device, dtype=torch.float64)
                    B = max(1, int(2e8 // (16 * eng.dim)))
                    for r0 in range(0, eng.dim, B):
                        blk = eng.psi[r0:r0 + B]
                        pa[r0:r0 + B] = (blk.real.double().square() + blk.imag.double().square()).sum(1)
                pa = pa.cpu().numpy()
                out = {"E_var": float(E), "p_string": pa / pa.sum(),
                       "strings": np.asarray(cistring.make_strings(range(norb), k), dtype=np.int64),
                       "t": time.time() - t0}
            elif req["cmd"] == "state_ci":
                # Task-1 protocol sampling (ffsim on the CPU state, C.task1_ci_strings); no SCI here
                import ffsim
                P = C.T1_PROTOCOL
                E = eng.energy(U, Z, t1=t1)
                v = eng.psi.detach().to("cpu").numpy().astype(np.complex128).reshape(-1)
                v /= np.linalg.norm(v)
                p = np.abs(v) ** 2
                ent = float(-(p[p > 0] * np.log(p[p > 0])).sum())
                p_hf = float(p[0])
                del p
                rng = np.random.default_rng(P["seed"])
                samples = ffsim.sample_state_vector(v, norb=norb, nelec=ham["nelec"], shots=P["shots_raw"],
                                                    seed=rng, bitstring_type=ffsim.BitstringType.INT)
                del v
                ci, nu = C.task1_ci_strings(samples, norb, ham["nelec"], rng, seed=P["seed"])
                out = {"E_var": float(E), "entropy": ent, "p_hf": p_hf, "n_unique_100k": nu,
                       "n_unique_1M": int(len(np.unique(samples))), "ci": ci, "t": time.time() - t0}
            else:
                out = {"error": f"unknown cmd {req.get('cmd')}"}
        except Exception as e:  # noqa: BLE001
            out = {"error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()}
        _write(proto_out, out)
    if sci is not None:
        sci.close()
    eng.release()


if __name__ == "__main__":
    main()
