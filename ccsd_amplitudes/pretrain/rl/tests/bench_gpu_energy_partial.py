#!/usr/bin/env python3
"""Component timing of pretrain/rl/gpu_energy.py for states that do not fit on this GPU (norb 18 on a
12 GB card), measured on partial states and scaled by exact counts:

  orbital rotation : minor-index Givens segments on a block of r rows (rows are independent) x dim/r,
                     4 minor rotations per energy (2 orbital rotations x (row, transposed row)),
  transpose / symmetrize : tiled square-matrix kernels on an m x m matrix, scaled by (dim/m)^2,
  phase passes     : the fp64 phase GEMM + polar + multiply on R rows x dim/R (2 per energy),
  energy           : one off-diagonal tile pair (2 gathers + 2 batched GEMMs + epilogue) and one diagonal
                     tile, from two (dim x T) column blocks, x the number of tile pairs.
Also run on norb 16/17 to check the estimate against the measured full-state time.

  CUDA_VISIBLE_DEVICES=2 python -m pretrain.rl.tests.bench_gpu_energy_partial --names C4H2_rxn9388_P \
      --tile 1408 --rows 2048 --out pretrain/rl/tests/results/bench_partial_n18.json
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import torch  # noqa: E402

from pretrain.rl.gpu_energy import LUCJEnergyGPU, givens_clustered, split_segments  # noqa: E402
from pretrain.rl.tests.bench_gpu_energy import haar  # noqa: E402


def timeit(f, rep=3):
    f()
    torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(rep):
        f()
    torch.cuda.synchronize()
    return (time.perf_counter() - t0) / rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--tile", type=int, default=None, help="energy tile T (default: what a 40 GB card would use)")
    ap.add_argument("--rows", type=int, default=2048)
    ap.add_argument("--square", type=int, default=16384)
    ap.add_argument("--target-mem-gb", type=float, default=40.0,
                    help="memory of the target card, used to pick the tile size it would use")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    res = []
    rng = np.random.default_rng(0)
    for name in args.names:
        torch.cuda.empty_cache()
        eng0 = LUCJEnergyGPU.from_npz(name, alloc_state=False, tile=64)   # tables only, to size things
        dim, K, NLp, n = eng0.dim, eng0.K, eng0.NLp, eng0.norb
        # tile a target card would choose: 55% of (budget - state - C table)
        itc = 8
        rest_t = args.target_mem_gb * 1e9 - 0.6e9 - itc * dim * dim - 4 * dim * NLp * K - 64e6
        T_target = int(math.sqrt(0.55 * max(rest_t, 1e8) / (itc * (K + 2 * NLp)))) // 32 * 32
        T = args.tile or min(T_target, 1024)                       # the engine's default max_tile
        eng0.release()
        del eng0
        torch.cuda.empty_cache()
        eng = LUCJEnergyGPU.from_npz(name, alloc_state=False, tile=T)
        be = eng.backend
        # ---- rotations: time on r1 and r2 rows, fit t(r) = a + b r (a = per-launch overhead), t(dim)
        W = haar(n, rng)
        app, D = givens_clustered(W, eng.t_top)
        segs = split_segments(app, n, eng.t_top)
        if eng.chunks is None:
            raise SystemExit("needs the numba fused (minor-index) backend")
        useg = be.upload_segments(segs)
        tr = {}
        for r in sorted({min(args.rows // 2, dim), min(args.rows * 2, dim)}):
            rows = torch.randn(r, dim, device="cuda", dtype=torch.complex64)
            rows_nb = be.arr(rows)
            tr[r] = timeit(lambda: be.givens_minor(rows_nb, r, useg, eng.pairs_b, eng.chunks))
            del rows, rows_nb
            torch.cuda.empty_cache()
        (r1, t1_), (r2, t2_) = sorted(tr.items())[0], sorted(tr.items())[-1]
        b_ = (t2_ - t1_) / max(r2 - r1, 1)
        a_ = max(t1_ - b_ * r1, 0.0)
        t_minor = a_ + b_ * dim                                             # one minor rotation, full state
        # ---- transpose / symmetrize
        m = min(args.square, dim)
        sq = torch.randn(m, m, device="cuda", dtype=torch.complex64)
        sq_nb = be.arr(sq)
        t_tr = timeit(lambda: be.transpose(sq_nb, m)) * (dim / m) ** 2
        t_sym = timeit(lambda: be.symmetrize(sq_nb, m)) * (dim / m) ** 2
        del sq, sq_nb
        torch.cuda.empty_cache()
        # ---- diagonal-Coulomb phase pass (numba kernel, square connectivity) on R rows
        R = min(dim, 2048)
        part = torch.randn(R, dim, device="cuda", dtype=torch.complex64)
        part_nb = be.arr(part)
        rowfac = torch.polar(torch.ones(dim, device="cuda", dtype=torch.float64),
                             torch.randn(dim, device="cuda", dtype=torch.float64))
        ez = torch.polar(torch.ones(n, device="cuda", dtype=torch.float64),
                         torch.randn(n, device="cuda", dtype=torch.float64))
        t_phase = timeit(lambda: be.phase_diag(part_nb, eng.g["strs_nb"], rowfac, ez, n, False)) * dim / R
        del part, part_nb
        torch.cuda.empty_cache()
        # ---- energy: one off-diagonal tile pair + one diagonal tile from two column blocks
        Tt = min(T, dim // 2)
        psiA = torch.randn(dim, Tt, device="cuda", dtype=torch.complex64)   # columns A = [0, Tt)
        psiB = torch.randn(dim, Tt, device="cuda", dtype=torch.complex64)   # columns B = [Tt, 2Tt)
        src = eng.g["src"]
        a0, a1, b0, b1 = 0, Tt, Tt, 2 * Tt

        def xt(cols, r0, r1, slot):
            nA, nB = r1 - r0, cols.shape[1]
            G = eng._Gbuf[: nA * K * nB].view(nA * K, nB)
            torch.index_select(cols, 0, src[r0:r1].reshape(-1), out=G)
            Gr = torch.view_as_real(G).view(nA, K, 2 * nB)
            off = slot * eng.T * NLp * eng.T
            Xc = eng._Xbuf[off: off + nA * NLp * nB]
            torch.bmm(eng.C[r0:r1], Gr, out=torch.view_as_real(Xc).view(nA, NLp, 2 * nB))
            return Xc.view(nA, NLp, nB)

        def pair():
            X1 = xt(psiB, a0, a1, 0)
            X2 = xt(psiA, b0, b1, 1)
            nb = be.epilogue(X1, X2, psiB[a0:a1], psiA[b0:b1], eng.nq, eng.npos, False, eng._out)
            eng._out[:nb].sum(0)

        def diag():
            X1 = xt(psiA, a0, a1, 0)
            nb = be.epilogue(X1, X1, psiA[a0:a1], psiA[a0:a1], eng.nq, eng.npos, True, eng._out)
            eng._out[:nb].sum(0)
        t_pair = timeit(pair)
        t_diag = timeit(diag)
        nt = -(-dim // T)
        n_off, n_diag = nt * (nt - 1) // 2, nt
        # every X tile (a-block, b-block) is computed once: cost per unit area of an off-diagonal pair x dim^2
        t_energy = (t_pair / (2 * Tt * Tt)) * dim * dim if nt > 1 else t_diag
        del psiA, psiB
        est = dict(rot=4 * t_minor + 2 * t_tr, phases=2 * t_phase, symmetrize=t_sym, energy=t_energy)
        est["total"] = sum(est.values())
        r_ = dict(name=name, norb=n, k=eng.k, dim=dim, t_top=eng.t_top, R_fused=eng.chunks["R"], K=K, NLp=NLp,
                  tile=T, tile_target=T_target, minor_rotation_s=t_minor, rot_rows_timing={str(k): v for k, v in tr.items()},
                  rot_launch_overhead_s=a_, transpose_s=t_tr, symmetrize_s=t_sym,
                  phase_pass_s=t_phase, tile_pair_s=t_pair, tile_diag_s=t_diag, n_tile_pairs=n_off,
                  n_diag_tiles=n_diag, estimate=est, mem=dict(psi_gb=8 * dim * dim / 1e9,
                                                                 ctab_gb=4 * dim * NLp * K / 1e9),
                  gpu=torch.cuda.get_device_name())
        res.append(r_)
        print(f"{name} n={n} k={eng.k} dim={dim} t={eng.t_top} T={T}: minor rot {t_minor:.2f}s, transpose "
              f"{t_tr:.3f}s, phase {t_phase:.3f}s, tile pair {1e3*t_pair:.1f} ms x {n_off} + diag "
              f"{1e3*t_diag:.1f} ms x {n_diag} -> energy {t_energy:.1f}s | estimate "
              + " ".join(f"{k}={v:.2f}" for k, v in est.items())
              + f" | psi {8*dim*dim/1e9:.1f} GB, C {4*dim*NLp*K/1e9:.2f} GB", flush=True)
        eng.release()
        del eng
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
