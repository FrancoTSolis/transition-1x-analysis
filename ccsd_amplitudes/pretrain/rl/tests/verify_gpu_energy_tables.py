#!/usr/bin/env python3
"""[verifier] CPU-only structural checks of the fused-kernel tables and the Givens factorization of
pretrain/rl/gpu_energy.py at the production sizes (norb 15..18, every closed-shell shape in rhf_hamiltonians),
where no statevector reference exists.

  * chunk tables (plan_fusion's t for complex64 and complex128): chunks are contiguous, partition [0, dim),
    have one top-orbital pattern each, the right k_low, and their local uint16 pair tables mapped to global
    indices reproduce EXACTLY the global adjacent-pair table of every low rotation p (no missing / extra /
    swapped pairs, also at chunk boundaries); shared-memory request <= 48 KB.
  * givens_clustered(W, t): prod of the embedded 2x2 blocks (application order) times diag(D) == W, every
    block is SU(2), D unimodular, 'low' segments never touch the top t orbitals.
  * the string-space Givens action (x' = M x on (lo, hi) pairs) composes to Lambda^k(W) on a small space,
    compared with the k x k minors of W (exact definition of the orbital rotation on Slater determinants).
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from pretrain.rl.gpu_energy import (chunk_tables, givens_clustered, plan_fusion, split_segments,  # noqa: E402
                                    string_tables)


def haar(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def check_chunks(norb, k, t, itemsize):
    tb = string_tables(norb, k)
    ct = chunk_tables(norb, k, t)
    strs = tb["strs"]
    dim = tb["dim"]
    nl = norb - t
    starts, lens, kls = ct["starts"], ct["lens"], ct["kls"]
    assert starts[0] == 0 and int(starts[-1] + lens[-1]) == dim
    assert np.all(starts[1:] == starts[:-1] + lens[:-1])
    for s, L, kl in zip(starts, lens, kls):
        top = strs[s:s + L] >> nl
        assert np.all(top == top[0])
        assert kl == k - bin(int(top[0])).count("1")
        assert L == comb(nl, int(kl))
        low = strs[s:s + L] & ((1 << nl) - 1)
        assert np.all(np.diff(low) > 0)
    n_pairs_checked = 0
    for p in range(nl - 1):                       # 'low' rotations: p + 1 < norb - t
        glo, ghi = tb["pairs"][p]
        gset = set(zip(glo.tolist(), ghi.tolist()))
        mine = set()
        for s, kl in zip(starts, kls):
            o, c = int(ct["off"][kl, p]), int(ct["cnt"][kl, p])
            lo16 = ct["lo16"][o:o + c].astype(np.int64)
            hi16 = ct["hi16"][o:o + c].astype(np.int64)
            # the same through the kernel's relative offsets into the per-k_low block
            o2 = int(ct["kl_base"][kl] + ct["off_rel"][kl, p])
            assert o2 == o
            assert np.array_equal(lo16, ct["lo"][o:o + c]) and np.array_equal(hi16, ct["hi"][o:o + c])
            for a, b in zip((s + lo16).tolist(), (s + hi16).tolist()):
                assert (a, b) not in mine
                mine.add((a, b))
        assert mine == gset, (norb, k, t, p, len(mine), len(gset))
        n_pairs_checked += len(gset)
    stride = int(ct["max_len"])
    pair_b = 4 * ct["max_kl_len"] + 16
    R = max(1, min(16, (48 * 1024 - pair_b) // (stride * itemsize)))
    smem = R * stride * itemsize + pair_b
    return dict(norb=norb, k=k, t=t, itemsize=itemsize, nch=len(starts), max_len=stride, R=int(R),
                smem_bytes=int(smem), smem_ok=bool(smem <= 48 * 1024), low_pairs_checked=n_pairs_checked)


def check_givens(n, t, rng, reps=3):
    worst = 0.0
    for _ in range(reps):
        W = haar(n, rng)
        app, D = givens_clustered(W, t)
        M = np.eye(n, dtype=complex)
        for p, B in app:               # application order: D first, then app[0], app[1], ...
            G = np.eye(n, dtype=complex)
            G[p:p + 2, p:p + 2] = B
            M = G @ M
            assert abs(np.linalg.det(B) - 1) < 1e-12 and np.abs(B @ B.conj().T - np.eye(2)).max() < 1e-12
        M = M @ np.diag(D)
        worst = max(worst, float(np.abs(M - W).max()))
        assert np.abs(np.abs(D) - 1).max() < 1e-12
        for kind, rots in split_segments(app, n, t):
            for p, _ in rots:
                assert (p + 1 < n - t) == (kind == "low")
    return worst


def lam_k(W, norb, k):
    """Lambda^k(W)[S, T] = det W[S, T] (rows = new string S, cols = old string T) in pyscf string order."""
    from pyscf.fci import cistring
    strs = cistring.make_strings(range(norb), k)
    occ = [[p for p in range(norb) if (s >> p) & 1] for s in strs]
    dim = len(strs)
    A = np.zeros((dim, dim), dtype=complex)
    for i, S in enumerate(occ):
        for j, T in enumerate(occ):
            A[i, j] = np.linalg.det(W[np.ix_(S, T)])
    return A


def check_string_action(norb, k, t, rng):
    tb = string_tables(norb, k)
    W = haar(norb, rng)
    app, D = givens_clustered(W, t)
    dim = tb["dim"]
    A = np.eye(dim, dtype=complex)
    # D first: string factor prod_{p in S} D_p
    A = np.diag(np.prod(np.where(tb["occ"], D[None, :], 1.0), axis=1)) @ A
    for p, M in app:
        lo, hi = tb["pairs"][p]
        x, y = A[lo].copy(), A[hi].copy()
        A[lo] = M[0, 0] * x + M[0, 1] * y
        A[hi] = M[1, 0] * x + M[1, 1] * y
    return float(np.abs(A - lam_k(W, norb, k)).max())


def main():
    rng = np.random.default_rng(99)
    out = dict(chunks=[], givens=[], string_action=[])
    shapes = [(15, 8), (15, 7), (16, 8), (16, 9), (16, 7), (17, 9), (17, 10), (17, 8), (18, 9), (18, 10),
              (18, 11), (18, 8), (12, 6), (13, 5), (14, 8)]
    for norb, k in shapes:
        for itemsize in (8, 16):
            plan = plan_fusion(norb, k, itemsize)
            if plan is None:
                out["chunks"].append(dict(norb=norb, k=k, itemsize=itemsize, plan=None))
                print(f"({norb},{k}) itemsize {itemsize}: no fusion plan", flush=True)
                continue
            t, R = plan
            r = check_chunks(norb, k, t, itemsize)
            r["plan_R"] = R
            out["chunks"].append(r)
            print(f"({norb},{k}) itemsize {itemsize}: t={t} R={R} nch={r['nch']} max_len={r['max_len']} "
                  f"smem={r['smem_bytes']} ok={r['smem_ok']}  low pairs checked {r['low_pairs_checked']}", flush=True)
        # an extra non-default split
        t2 = min(norb - 2, (plan_fusion(norb, k, 8) or (3, 0))[0] + 1)
        r = check_chunks(norb, k, t2, 8)
        out["chunks"].append(r)
        print(f"({norb},{k}) extra t={t2}: nch={r['nch']} ok={r['smem_ok']} pairs {r['low_pairs_checked']}", flush=True)
    for n in (15, 16, 17, 18):
        for t in (0, 1, 3, 5, 6, 7):
            w = check_givens(n, t, rng)
            out["givens"].append(dict(n=n, t=t, max_abs_err=w))
        print(f"givens_clustered n={n}: max |prod - W| = {max(g['max_abs_err'] for g in out['givens'] if g['n'] == n):.1e}",
              flush=True)
    for norb, k, t in [(6, 3, 2), (7, 2, 3), (7, 5, 1), (8, 4, 3), (8, 3, 5), (9, 4, 4), (9, 6, 2)]:
        e = check_string_action(norb, k, t, rng)
        out["string_action"].append(dict(norb=norb, k=k, t=t, max_abs_err=e))
        print(f"string action vs Lambda^k(W) minors: ({norb},{k}) t={t}: {e:.1e}", flush=True)
    p = ROOT / "pretrain/rl/tests/results/verify_gpu_energy_tables.json"
    json.dump(out, open(p, "w"), indent=1)
    print("all structural checks passed ->", p)


if __name__ == "__main__":
    main()
