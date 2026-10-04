#!/usr/bin/env python3
"""[verifier] Independent random validation of pretrain/rl/gpu_energy.py against ffsim.

New seeds and edge cases the builder's test did not cover: norb 2..14, k = 1 and k = norb-1, nocc > nvirt,
n_reps 1/2/3, t1 None / given, non-unitary U (polar), square / hex / heavy-hex / all-to-all, Hamiltonians =
random sym8 / random PSD / sub-blocks of the real rhf_hamiltonians integrals (realistic ERI spectra);
engine variants: default, forced small tiles (multi-tile path), forced fusion split t and tiny shared memory
(different R / chunking), fuse=False, torch backend.  Also energy_of_state on random complex symmetric states
(generic states, not LUCJ) vs ffsim's linear operator Rayleigh quotient.

  CUDA_VISIBLE_DEVICES=2 OMP_NUM_THREADS=4 python -m pretrain.rl.tests.verify_gpu_energy_random --seed 777 --out X.json
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

from pretrain.rl.energy import exact_energy, make_ucj_op  # noqa: E402
from pretrain.rl.gpu_energy import LUCJEnergyGPU  # noqa: E402


def sym8(t):
    perms = [(0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2),
             (2, 3, 0, 1), (3, 2, 0, 1), (2, 3, 1, 0), (3, 2, 1, 0)]
    return sum(t.transpose(p) for p in perms) / 8


def random_ham(norb, rng, kind):
    h = rng.normal(size=(norb, norb))
    h = 0.5 * (h + h.T)
    if kind == "sym8":
        eri = sym8(rng.normal(size=(norb,) * 4)) * 0.7
    else:
        L = rng.normal(size=(norb * norb // 2 + 1, norb, norb))
        L = 0.5 * (L + L.transpose(0, 2, 1))
        eri = np.einsum("lpq,lrs->pqrs", L, L) / (norb * norb)
    return h, eri, float(rng.normal() * 3)


NPZ_NAMES = ["C2H3N_rxn2857_P", "C2H4O_rxn0724_R", "C3H4_rxn2391_TS", "C2N2_rxn3923_R", "C3HN_rxn2586_P"]


def npz_block(name, norb, rng):
    """Sub-block of a real molecule's integrals on a random subset of orbitals (keeps 8-fold symmetry and a
    realistic, nearly singular ERI spectrum)."""
    d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
    n0 = int(d["norb"])
    k0 = int(d["nelec_a"])
    # keep the orbitals around the Fermi level so the HF reference is sensible
    lo = max(0, k0 - norb // 2 - int(rng.integers(0, 2)))
    S = np.arange(lo, min(n0, lo + norb))
    if len(S) < norb:
        S = np.arange(n0 - norb, n0)
    h = d["one_body"][np.ix_(S, S)]
    eri = d["two_body"][np.ix_(S, S, S, S)]
    return h, eri, float(d["constant"])


def haar(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def near_id(n, rng, s):
    a = s * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    a = a - a.conj().T
    from scipy.linalg import expm
    return expm(a)


def polar(U):
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def build_cases(rng, n_cases):
    shapes = [(2, 1), (3, 1), (3, 2), (4, 1), (4, 3), (5, 1), (5, 4), (6, 1), (6, 5), (7, 6), (8, 2), (8, 7),
              (9, 2), (9, 7), (10, 2), (10, 8), (11, 3), (11, 8), (12, 4), (12, 6), (12, 9), (12, 3),
              (13, 6), (13, 9), (14, 7), (14, 5), (11, 6), (10, 5)]
    conns = ["square", "hex", "heavy-hex", "all-to-all", "square", "square"]
    hams = ["sym8", "psd", "npz", "npz"]
    cases = []
    for c in range(n_cases):
        norb, k = shapes[c % len(shapes)]
        cases.append(dict(norb=norb, k=k, conn=conns[c % len(conns)], ham=hams[c % len(hams)],
                          n_reps=[2, 1, 3, 2][c % 4], t1=["rand", "none", "rand"][c % 3],
                          Ukind=["haar", "near", "nonunit", "haar", "near"][c % 5]))
    return cases


def variants(norb, dim):
    v = [("numba_c16", dict(dtype=torch.complex128, backend="numba")),
         ("numba_c8", dict(dtype=torch.complex64, backend="numba")),
         ("numba_c16_tile32", dict(dtype=torch.complex128, backend="numba", tile=32)),
         ("numba_c8_tile32", dict(dtype=torch.complex64, backend="numba", tile=32)),
         ("numba_c16_nofuse_tile64", dict(dtype=torch.complex128, backend="numba", fuse=False, tile=64)),
         ("torch_c16_tile32", dict(dtype=torch.complex128, backend="torch", tile=32))]
    if norb >= 5:
        for t in sorted({1, 2, norb - 3, norb - 2}):
            if 1 <= t <= norb - 2:
                v.append((f"numba_c16_t{t}_smem6", dict(dtype=torch.complex128, backend="numba", fuse_t=t,
                                                       smem_kb=6, tile=64)))
        v.append(("numba_c8_smem3", dict(dtype=torch.complex64, backend="numba", smem_kb=3, tile=96)))
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=777)
    ap.add_argument("--n-cases", type=int, default=56)
    ap.add_argument("--n-generic", type=int, default=8)
    ap.add_argument("--max-norb", type=int, default=14)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import ffsim

    rng = np.random.default_rng(args.seed)
    rows, worst = [], {}
    cases = [c for c in build_cases(rng, args.n_cases) if c["norb"] <= args.max_norb]
    for ci, cs in enumerate(cases):
        norb, k, conn = cs["norb"], cs["k"], cs["conn"]
        if cs["ham"] == "npz":
            h, eri, const = npz_block(NPZ_NAMES[ci % len(NPZ_NAMES)], norb, rng)
            hk = f"npz:{NPZ_NAMES[ci % len(NPZ_NAMES)]}"
        else:
            h, eri, const = random_ham(norb, rng, cs["ham"])
            hk = cs["ham"]
        R = cs["n_reps"]
        if cs["Ukind"] == "haar":
            U = np.stack([haar(norb, rng) for _ in range(R)])
        elif cs["Ukind"] == "near":
            U = np.stack([near_id(norb, rng, 0.3) for _ in range(R)])
        else:  # non-unitary (like float32 network outputs, but larger): engine must apply the same polar
            U = np.stack([haar(norb, rng) + 3e-3 * (rng.normal(size=(norb, norb)) + 1j * rng.normal(size=(norb, norb)))
                          for _ in range(R)])
        Z = rng.normal(size=(R, norb, norb)) * float(rng.choice([0.2, 1.0, 2.5]))
        Z = 0.5 * (Z + Z.transpose(0, 2, 1))
        t1 = None
        if cs["t1"] == "rand" and norb - k > 0:
            t1 = float(rng.choice([0.1, 0.4])) * rng.normal(size=(k, norb - k))
        ham = ffsim.MolecularHamiltonian(h, eri, const)
        Up = polar(U)
        t0 = time.time()
        op = make_ucj_op(Z, Up, conn, t1=t1)
        E_ref = exact_energy(ham, norb, (k, k), op)
        psi_ref = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, (k, k)), op, norb=norb, nelec=(k, k))
        t_ref = time.time() - t0
        from math import comb
        dim = comb(norb, k)
        row = dict(case=ci, **cs, hamkey=hk, dim=dim, E_ref=E_ref, t_ref=t_ref, res={})
        for tag, kw in variants(norb, dim):
            try:
                eng = LUCJEnergyGPU(h, eri, const, norb, (k, k), max_mem_gb=3.0, **kw)
                # mix numpy / torch inputs
                if ci % 2 == 1:
                    Uin = torch.as_tensor(U, device="cuda")
                    tin = None if t1 is None else torch.as_tensor(t1)
                else:
                    Uin, tin = U, t1
                E = eng.energy(Uin, Z, t1=tin, connectivity=conn)
                psi = eng.state(U, Z, t1=t1, connectivity=conn).cpu().numpy().reshape(-1)
                dpsi = float(np.abs(psi - psi_ref).max())
                info = dict(E=E, dE=E - E_ref, dpsi=dpsi, T=eng.T, t_top=eng.t_top,
                            R=(eng.chunks or {}).get("R"), backend=eng.backend.name)
                eng.release()
                del eng
            except Exception as e:  # noqa: BLE001
                info = dict(error=f"{type(e).__name__}: {e}")
            row["res"][tag] = info
            if "error" not in info:
                w = worst.setdefault(tag, dict(dE=0.0, dpsi=0.0, n=0))
                w["dE"] = max(w["dE"], abs(info["dE"]))
                w["dpsi"] = max(w["dpsi"], info["dpsi"])
                w["n"] += 1
        rows.append(row)
        errs = [t for t, r in row["res"].items() if "error" in r]
        c16 = max(abs(r["dE"]) for t, r in row["res"].items() if "error" not in r and "c16" in t)
        c8 = max(abs(r["dE"]) for t, r in row["res"].items() if "error" not in r and "c8" in t)
        dps = max(r["dpsi"] for t, r in row["res"].items() if "error" not in r and "c16" in t)
        print(f"[{ci:2d}] n={norb:2d} k={k:2d} dim={dim:5d} {conn:10s} {hk:22s} reps={R} t1={cs['t1']:4s} U={cs['Ukind']:7s}"
              f" E_ref {E_ref:+.10f}  max|dE| c16 {c16:.1e} c8 {c8:.1e}  max dpsi(c16) {dps:.1e}"
              f"  ({len(row['res'])} variants, ref {t_ref:.1f}s){'  ERR ' + str(errs) if errs else ''}", flush=True)
        for t in errs:
            print("      ", t, row["res"][t]["error"], flush=True)

    # generic symmetric complex states through energy_of_state (multi-tile), vs ffsim Rayleigh quotient
    gen = []
    gshapes = [(8, 4), (10, 3), (11, 5), (12, 7), (13, 6), (14, 7), (12, 5), (9, 6)][: args.n_generic]
    for gi, (norb, k) in enumerate(gshapes):
        if norb > args.max_norb:
            continue
        if gi % 2 == 0:
            h, eri, const = npz_block(NPZ_NAMES[gi % len(NPZ_NAMES)], norb, rng)
        else:
            h, eri, const = random_ham(norb, rng, "sym8")
        from math import comb
        dim = comb(norb, k)
        A = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        psi = (A + A.T) / 2
        # mix in a dominant HF-like component so it resembles physical states, half of the time
        if gi % 2 == 0:
            psi = 0.05 * psi / np.linalg.norm(psi)
            psi[0, 0] += 1.0
        psi = psi / np.linalg.norm(psi) * (1.0 + 0.1 * gi)          # unnormalized on purpose
        ham = ffsim.MolecularHamiltonian(h, eri, const)
        lin = ffsim.linear_operator(ham, norb=norb, nelec=(k, k))
        v = psi.reshape(-1)
        E_ref = float(np.vdot(v, lin @ v).real / np.vdot(v, v).real)
        r = dict(norb=norb, k=k, dim=dim, E_ref=E_ref, res={})
        for tag, kw in [("c16_tile32", dict(dtype=torch.complex128, tile=32)),
                        ("c16_default", dict(dtype=torch.complex128)),
                        ("c8_tile64", dict(dtype=torch.complex64, tile=64)),
                        ("torch_c16_tile32", dict(dtype=torch.complex128, backend="torch", tile=32))]:
            eng = LUCJEnergyGPU(h, eri, const, norb, (k, k), max_mem_gb=3.0, **kw)
            E = eng.energy_of_state(psi)
            r["res"][tag] = dict(E=E, dE=E - E_ref, T=eng.T)
            eng.release()
        gen.append(r)
        print(f"[gen {gi}] n={norb} k={k} dim={dim} E_ref {E_ref:+.10f}  " +
              "  ".join(f"{t}: dE {x['dE']:+.1e} (T={x['T']})" for t, x in r["res"].items()), flush=True)

    print("\n=== worst over", len(rows), "LUCJ cases ===")
    for t, w in worst.items():
        print(f"  {t:28s} n={w['n']:3d}  max|dE| {w['dE']:.2e}  max|dpsi| {w['dpsi']:.2e}")
    if gen:
        print("=== generic states ===")
        for t in gen[0]["res"]:
            print(f"  {t:20s} max|dE| {max(abs(g['res'][t]['dE']) for g in gen):.2e}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(worst=worst, rows=rows, generic=gen, args=vars(args)), open(args.out, "w"), indent=1)
    print("->", args.out)


if __name__ == "__main__":
    main()
