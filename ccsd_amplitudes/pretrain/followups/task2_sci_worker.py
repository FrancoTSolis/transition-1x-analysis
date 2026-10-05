#!/usr/bin/env python3
"""Selected-CI (QSCI) server for task 2: a pure-CPU process that never imports torch.

Why a separate process: in pretrain/.train_venv, PySCF's selected CI returned wrong subspace energies (both the
Davidson eigenvalue and the RDM energy, e.g. -130.31293 / -130.32388 instead of -130.38340 Ha) in processes where
`pyscf.lib` had been imported before torch and the GPU engine (two vendored libgomp copies are loaded; the exact
symbol mix-up was not pinned down).  The values from torch-free processes are thread-count independent and agree
to 1e-12 Ha with the full-space Rayleigh quotient of the returned eigenvector on the GPU engine (complex128), and
the GPU worker re-checks every objective value that way (complex64).

Requests (length-prefixed pickles on stdin/stdout, see task2_worker):
  {"cmd": "qsci", "bitstrings", "settings", "return_vec"} -> task2_common.qsci_energy (+ best eigenvector)
      the subsampling rng is persistent across requests (default_rng(entropy)), as in Lin et al.'s optimization
  {"cmd": "solve", "ci": [(strs_a, strs_b), ...], "spin_sq"} -> energies of the given subspaces (identical
      subspaces solved once), + the eigenvector of the lowest one when return_vec
  {"cmd": "quit"}
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time
import traceback
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--entropy", type=int, default=0)
    args = ap.parse_args()
    proto_out = os.fdopen(os.dup(1), "wb")
    os.dup2(2, 1)
    proto_in = sys.stdin.buffer
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import numpy as np
    from pyscf import lib
    from qiskit_addon_sqd.fermion import solve_sci
    lib.num_threads(args.threads)
    from pretrain.followups import task2_common as C
    from pretrain.followups.task2_worker import _read, _write
    assert "torch" not in sys.modules, "the SCI worker must not import torch"
    ham = C.load_ham(args.name)
    rng = np.random.default_rng(args.entropy)
    _write(proto_out, {"ready": True, "pid": os.getpid()})
    while True:
        req = _read(proto_in)
        if req is None or req.get("cmd") == "quit":
            break
        t0 = time.time()
        try:
            if req["cmd"] == "qsci":
                q = C.qsci_energy(ham, req["bitstrings"], req["settings"], rng, keep_best=req.get("return_vec", False))
                best = q.pop("best", None)
                if best is not None:
                    q["vec"] = (np.asarray(best.sci_state.amplitudes), np.asarray(best.sci_state.ci_strs_a),
                                np.asarray(best.sci_state.ci_strs_b))
                out = q
            elif req["cmd"] == "solve":
                memo, Es, best = {}, [], None
                for sa, sb in req["ci"]:
                    key = hashlib.sha1(np.ascontiguousarray(sa).tobytes() + b"|" +
                                       np.ascontiguousarray(sb).tobytes()).hexdigest()
                    if key not in memo:
                        r = solve_sci((sa, sb), ham["one_body"], ham["two_body"], norb=ham["norb"],
                                      nelec=ham["nelec"], spin_sq=req.get("spin_sq"))
                        e = float(r.energy + ham["constant"])
                        memo[key] = e
                        if best is None or e < best[0]:
                            best = (e, r)
                    Es.append(memo[key])
                out = {"E_batches": Es, "n_distinct": len(memo)}
                if req.get("return_vec") and best is not None:
                    st = best[1].sci_state
                    out["vec"] = (np.asarray(st.amplitudes), np.asarray(st.ci_strs_a), np.asarray(st.ci_strs_b))
                    out["E_vec"] = best[0]
            else:
                out = {"error": f"unknown cmd {req.get('cmd')}"}
            out["t_sci"] = time.time() - t0
        except Exception as e:  # noqa: BLE001
            out = {"error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()}
        _write(proto_out, out)


if __name__ == "__main__":
    main()
