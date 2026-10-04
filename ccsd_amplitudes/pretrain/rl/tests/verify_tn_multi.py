#!/usr/bin/env python3
"""[Independent verifier] several LUCJEnergyTN evaluators (different molecules) alive in ONE process.

block2's DMRGDriver replaces the global data frame (Global.frame) on construction.  This checks whether an
evaluator still returns correct energies after another evaluator was created / deleted (typical RL worker: one
evaluator per molecule), including a cache miss (new t1 -> MPO built while another driver's frame is global).
Usage: python3 pretrain/rl/tests/verify_tn_multi.py
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
    os.environ[_v] = "2"
import gc  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import ffsim  # noqa: E402

from pretrain.rl import tn_energy as T  # noqa: E402
from pretrain.rl.energy import exact_energy, make_ucj_op  # noqa: E402

SCR = os.environ.get("TN_VERIFY_SCRATCH", tempfile.gettempdir())


def random_unitary(n, rng):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(X)
    return Q * (np.diagonal(R) / abs(np.diagonal(R)))


def system(n, nelec, seed):
    rng = np.random.default_rng(seed)
    ham = ffsim.random.random_molecular_hamiltonian(n, seed=seed, dtype=float)
    U = np.stack([random_unitary(n, rng) for _ in range(2)])
    Z = rng.normal(scale=0.5, size=(2, n, n))
    Z = 0.5 * (Z + Z.transpose(0, 2, 1))
    t1a = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
    t1b = rng.normal(scale=0.1, size=(nelec[0], n - nelec[0]))
    ref = {tag: exact_energy(ham, n, nelec, make_ucj_op(Z, U, "square", t1=t)) for tag, t in (("a", t1a), ("b", t1b))}
    return ham, U, Z, {"a": t1a, "b": t1b}, ref


def make_ev(ham, n, nelec):
    return T.LUCJEnergyTN(ham.one_body_tensor, ham.two_body_tensor, ham.constant, n, nelec, max_bond=10 ** 6,
                          device="cpu", basis="er", cutoff=0.0, block2_threads=1, stack_mem=1 << 28,
                          scratch=tempfile.mkdtemp(prefix="vmulti_", dir=SCR))


def check(tag, ev, sysd, t1key):
    ham, U, Z, t1s, ref = sysd
    E, _ = ev.energy(U, Z, t1s[t1key])
    d = E - ref[t1key]
    print(f"{tag:60s} E {E:.10f} ref {ref[t1key]:.10f} diff {d:+.2e}  {'OK' if abs(d) < 1e-7 else 'MISMATCH'}",
          flush=True)
    return d


if __name__ == "__main__":
    A = system(6, (3, 3), 101)
    B = system(7, (4, 4), 202)
    C = system(5, (2, 2), 303)
    diffs = []
    evA = make_ev(A[0], 6, (3, 3))
    diffs.append(check("1. evA alone, t1 a", evA, A, "a"))
    evB = make_ev(B[0], 7, (4, 4))
    diffs.append(check("2. evB created, evB t1 a", evB, B, "a"))
    diffs.append(check("3. evA again (cached MPO from own frame), t1 a", evA, A, "a"))
    diffs.append(check("4. evA new t1 b (MPO built while evB frame is global)", evA, A, "b"))
    del evB
    gc.collect()
    evC = make_ev(C[0], 5, (2, 2))
    diffs.append(check("5. evB deleted, evC created, evC t1 a", evC, C, "a"))
    diffs.append(check("6. evA t1 b (cached MPO allocated in evB's freed frame?)", evA, A, "b"))
    diffs.append(check("7. evA t1 a (rebuild)", evA, A, "a"))
    evB2 = make_ev(B[0], 7, (4, 4))
    diffs.append(check("8. evB2 (re-created) t1 b", evB2, B, "b"))
    diffs.append(check("9. evC t1 b", evC, C, "b"))
    diffs.append(check("10. evA t1 a", evA, A, "a"))
    print("max |diff|", max(abs(d) for d in diffs))
