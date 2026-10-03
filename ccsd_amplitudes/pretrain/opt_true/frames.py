#!/usr/bin/env python3
"""Chemistry frames for the masked compressed-DF (optimize=True LUCJ) problem.

Finding behind this module (design panel + own checks): the optimize=True orbitals are atom-centred hybrids
(occupied weight ~0.5, 1.4-2.2 effective atoms), and every optimized column has a relative occupied/virtual phase
of -pi/4 (rep 0) or +pi/4 (rep 1).  So we work in the REAL SECTOR
    U_k = Phi_k B_k expm(A_k),   Phi_0 = diag(1_occ, e^{-i pi/4} 1_vir),  Phi_1 = conj(Phi_0),  A_k real antisym.
where B_k is a deterministic, geometry-defined orthogonal frame of the active space (STO-3G valence = active space):
  * OAOs: valence AOs (heavy-atom 1s excluded) projected on the active MOs, Loewdin-orthonormalized
    (MO-sign covariant: B rows flip with the MO signs).
  * Canonical orientation: p orbitals are taken along the molecule's principal axes (sign-fixed by third moments,
    right-handed), so the frame does not depend on the lab orientation of the input geometry.
  * Rep 0: bond-directed sp3-like hybrids on heavy atoms (lone pairs placed at tetrahedral angles from bond
    geometry); rep 1: plain OAOs (s, p1, p2, p3 along the principal axes).  This breaks the rep-swap symmetry.
  * Chain order (the square mask couples p, p+1): DFS over heavy atoms with every tree bond made chain-adjacent:
    [A->parent, (A->H_i, H_i s)..., lone pairs of A, (A->child_j, subtree of child_j)...].  Rep 1 assigns the
    OAOs of each atom to the same slots (slot -> atom map identical across reps).

Output (precompute): rhf_frames/<nocc>_<nvirt>.npz with names, B (M,2,n,n) float32, slot_atom (M,n) int16,
slot_type (M,2,n) int8, atom_Z (M,natm_max) int8, coords (M,natm_max,3) float32 (canonical orientation, Angstrom),
bond (M,natm_max,natm_max) bool, natm (M,).

Usage: python3 -m pretrain.opt_true.frames --n-procs 24 [--shapes 16_13]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57}
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998}
ZNUM = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9}
# slot types: 0 H-1s | rep0: 1 hybrid->parent, 2 hybrid->H, 3 lone pair, 4 hybrid->child, 5 isolated-atom AO
#                    | rep1: 6 s, 7 p1, 8 p2, 9 p3
N_SLOT_TYPES = 10
TET = np.arccos(-1.0 / 3.0)


def read_xyz(name):
    L = open(ROOT / "jobs" / name / f"{name}.xyz").read().splitlines()
    k = int(L[0])
    sym, R = [], []
    for ln in L[2:2 + k]:
        p = ln.split()
        sym.append(p[0])
        R.append([float(x) for x in p[1:4]])
    return sym, np.array(R)


def principal_axes(sym, R):
    """Columns = canonical axes (lab coordinates): largest spread first, signs by third moments, right-handed."""
    m = np.array([MASS[s] for s in sym])
    X = R - (m[:, None] * R).sum(0) / m.sum()
    C2 = np.einsum("i,ij,ik->jk", m, X, X)
    w, V = np.linalg.eigh(C2)
    V = V[:, ::-1].copy()
    for j in range(2):
        proj = X @ V[:, j]
        t = (m * proj ** 3).sum()
        if abs(t) < 1e-6:                       # symmetric along this axis: fall back to the heaviest-atom side
            t = proj[np.argmax(m + 1e-3 * np.arange(len(m)))]
        if t < 0:
            V[:, j] *= -1
    V[:, 2] = np.cross(V[:, 0], V[:, 1])
    return V


def bonds(sym, R, scale=1.25):
    rad = np.array([COV[s] for s in sym])
    D = np.linalg.norm(R[:, None] - R[None], axis=-1)
    return (D < scale * (rad[:, None] + rad[None])) & (D > 0), D


def _hyb(d):
    return np.r_[1.0, np.sqrt(3.0) * d] / 2.0


def _perp(d, refs):
    for r in refs:
        u = r - (r @ d) * d
        if np.linalg.norm(u) > 1e-3:
            return u / np.linalg.norm(u)
    return None


def _lowdin(V):
    w, Q = np.linalg.eigh(V.T @ V)
    return V @ Q @ np.diag(w ** -0.5) @ Q.T, w.min()


def atom_hybrids(A, nbs, Rc, axes_c, heavy_nb):
    """4 orthonormal hybrid vectors (s,p1,p2,p3 coefficients, canonical axes) of heavy atom A.

    nbs: bonded partners in priority order (parent, H atoms, other heavy atoms).
    1. bond hybrids (s + sqrt3 d.p)/2 toward up to 4 partners, Loewdin-orthonormalized among themselves
       (partners are dropped, farthest first, while the set is ill-conditioned: TS / near-linear geometry);
    2. lone pairs fill the exact orthogonal complement: geometric candidates (tetrahedral placement from the bond
       directions, then the pure p orbital normal to the bond plane, then s, p1, p2, p3) are projected on the
       complement and accepted by Gram-Schmidt.  Planar 3-bond atoms thus get sp2-like bond hybrids + p_pi.
    Returns list of (partner or None, vec); None = lone pair."""
    unit = lambda B: (Rc[B] - Rc[A]) / np.linalg.norm(Rc[B] - Rc[A])  # noqa: E731
    nbs = list(nbs)
    if len(nbs) > 4:
        keep = set(sorted(nbs, key=lambda B: np.linalg.norm(Rc[B] - Rc[A]))[:4])
        nbs = [B for B in nbs if B in keep]
    while nbs:
        Vb, wmin = _lowdin(np.array([_hyb(unit(B)) for B in nbs]).T)
        if wmin > 0.05:
            break
        far = max(nbs, key=lambda B: np.linalg.norm(Rc[B] - Rc[A]))
        nbs = [B for B in nbs if B != far]
    if not nbs:
        return [(None, e) for e in np.eye(4)]
    dirs = [unit(B) for B in nbs]
    k = len(dirs)
    cands = []
    if k == 3:
        sd = sum(dirs)
        if np.linalg.norm(sd) > 1e-3:
            cands.append(_hyb(-sd / np.linalg.norm(sd)))
        nrm = np.cross(dirs[0] - dirs[2], dirs[1] - dirs[2])
        if np.linalg.norm(nrm) > 1e-6:
            nrm = nrm / np.linalg.norm(nrm)
            cands.append(np.r_[0.0, nrm if nrm @ axes_c[:, 2] >= 0 else -nrm])
    elif k == 2:
        bsec = -(dirs[0] + dirs[1])
        bsec = bsec / np.linalg.norm(bsec) if np.linalg.norm(bsec) > 1e-3 else _perp(dirs[0], list(axes_c.T))
        m = np.cross(dirs[0], dirs[1])
        m = _perp(bsec, [m] + list(axes_c.T)) if np.linalg.norm(m) > 1e-3 else _perp(bsec, list(axes_c.T))
        if m @ axes_c[:, 2] < 0 or (abs(m @ axes_c[:, 2]) < 1e-3 and m @ axes_c[:, 1] < 0):
            m = -m
        c, s_ = np.cos(TET / 2), np.sin(TET / 2)
        cands += [_hyb(c * bsec + s_ * m), _hyb(c * bsec - s_ * m)]
    elif k == 1:
        d = dirs[0]
        refs = [unit(B) for B in heavy_nb.get(nbs[0], []) if B != A] + list(axes_c.T)   # dihedral reference first
        u = _perp(d, refs)
        w = np.cross(d, u)
        cands += [_hyb(-d / 3.0 + (2 * np.sqrt(2) / 3.0) * (np.cos(ph) * u + np.sin(ph) * w))
                  for ph in (0.0, 2 * np.pi / 3, 4 * np.pi / 3)]
    cands += list(np.eye(4))
    P = np.eye(4) - Vb @ Vb.T
    lps = []
    for c in cands:
        if len(lps) == 4 - k:
            break
        v = P @ c
        for l in lps:
            v = v - (l @ v) * l
        if np.linalg.norm(v) > 0.3 * np.linalg.norm(c):
            lps.append(v / np.linalg.norm(v))
    if len(lps) > 1 and k >= 1:                       # symmetric orthonormalization within the complement
        L, _ = _lowdin(np.array([P @ c for c in cands[:len(lps)]]).T) if all(
            np.linalg.norm(P @ c) > 0.3 * np.linalg.norm(c) for c in cands[:len(lps)]) else (np.array(lps).T, 1.0)
        if np.abs(L.T @ L - np.eye(len(lps))).max() < 1e-8 and np.abs(Vb.T @ L).max() < 1e-8:
            lps = list(L.T)
    parts = list(nbs) + [None] * len(lps)
    V = np.concatenate([Vb, np.array(lps).T], 1) if lps else Vb
    assert V.shape == (4, 4) and np.abs(V.T @ V - np.eye(4)).max() < 1e-8, (A, k, len(lps))
    return [(pt, V[:, j]) for j, pt in enumerate(parts)]


def mayer_bonds(mol, C, nocc, sym, thresh):
    """Mayer bond orders from the valence part of the RHF density (active occupied MOs)."""
    S = mol.intor("int1e_ovlp")
    Cocc = C[:, :nocc]
    PS = 2.0 * Cocc @ Cocc.T @ S
    labs = mol.ao_labels(fmt=False)
    ao_atom = np.array([l[0] for l in labs])
    nat = len(sym)
    BO = np.zeros((nat, nat))
    for A in range(nat):
        ia = ao_atom == A
        for B in range(A + 1, nat):
            ib = ao_atom == B
            BO[A, B] = BO[B, A] = float((PS[np.ix_(ia, ib)] * PS[np.ix_(ib, ia)].T).sum())
    return BO > thresh, BO


def build_frame(name, bond_scale=1.25, bond_mode="distance", bo_thresh=0.3):
    from pyscf import gto
    sym, R = read_xyz(name)
    mol = gto.M(atom=[(s, tuple(r)) for s, r in zip(sym, R)], basis="sto-3g", unit="Angstrom", verbose=0)
    z = np.load(ROOT / "rhf_dataset" / f"{name}.npz")
    C = z["mo_coeff"].astype(np.float64)
    S = mol.intor("int1e_ovlp")
    labs = mol.ao_labels(fmt=False)
    val = [mu for mu, l in enumerate(labs) if not (sym[l[0]] != "H" and l[2] == "1s")]
    n = C.shape[1]
    assert len(val) == n, (name, len(val), n)
    M = (C.T @ S)[:, val]
    w, V = np.linalg.eigh(M.T @ M)
    O = M @ V @ np.diag(w ** -0.5) @ V.T                                     # (n MO, n OAO), orthogonal
    ax = principal_axes(sym, R)
    Rc = (R - R.mean(0)) @ ax                                                 # canonical coordinates
    axes_c = np.eye(3)
    # per-atom OAO columns in canonical orientation: (s, p1, p2, p3) or (s,)
    ao_atom = np.array([labs[mu][0] for mu in val])
    cols = {}
    for A in range(len(sym)):
        idx = np.where(ao_atom == A)[0]
        if sym[A] == "H":
            cols[A] = O[:, idx]
        else:
            s_i = [i for i in idx if labs[val[i]][2].endswith("s")]
            p_i = sorted([i for i in idx if labs[val[i]][2].endswith("p")], key=lambda i: "xyz".index(labs[val[i]][3]))
            p_lab = O[:, p_i]                                                 # px, py, pz
            cols[A] = np.concatenate([O[:, s_i], p_lab @ ax], 1)              # s, p1, p2, p3 (canonical)
    bond, D = bonds(sym, R, bond_scale)
    if bond_mode in ("mayer", "either"):
        nocc = int(z["nocc"]) if "nocc" in z.files else None
        mb, _ = mayer_bonds(mol, C, nocc, sym, bo_thresh)
        bond = mb if bond_mode == "mayer" else (bond | mb)
    heavy = [A for A in range(len(sym)) if sym[A] != "H"]
    heavy_nb = {A: [B for B in heavy if bond[A, B]] for A in heavy}
    # H atoms: attach each H to its closest bonded heavy atom (or closest heavy atom if none bonded)
    h_owner = {}
    for Hh in [A for A in range(len(sym)) if sym[A] == "H"]:
        cand = [A for A in heavy if bond[Hh, A]] or heavy
        h_owner[Hh] = min(cand, key=lambda A: D[Hh, A])
    h_of = defaultdict(list)
    for Hh, A in h_owner.items():
        h_of[A].append(Hh)
    for A in h_of:
        h_of[A].sort(key=lambda Hh: -Rc[Hh, 0])
    # DFS start: terminal heavy atom (fewest heavy neighbours), largest coordinate along the long axis
    seen, slots0, slot_atom, slots1_assign = set(), [], [], {}

    def emit(kind, partner, vec, A):
        slots0.append((kind, partner, cols[A] @ vec))
        slot_atom.append(A)

    def emit_h(Hh):
        slots0.append(("Hs", Hh, cols[Hh][:, 0]))
        slot_atom.append(Hh)

    def dfs(A, parent):
        seen.add(A)
        others = sorted([B for B in heavy_nb[A] if B != parent], key=lambda B: (B in seen, -Rc[B, 0]))
        nbs = ([parent] if parent is not None else []) + list(h_of[A]) + others
        hyb = atom_hybrids(A, nbs, Rc, axes_c, heavy_nb)
        by_partner = {pt: v for pt, v in hyb if pt is not None}
        start = len(slot_atom)
        if parent is not None and parent in by_partner:
            emit("parent", parent, by_partner[parent], A)
        for Hh in h_of[A]:
            if Hh in by_partner:
                emit("H", Hh, by_partner[Hh], A)
            emit_h(Hh)
        for pt, v in hyb:
            if pt is None:
                emit("lp", None, v, A)
        for B in others:
            if B in by_partner:
                emit("child", B, by_partner[B], A)
            if B not in seen:
                dfs(B, A)
        own = [i for i in range(start, len(slot_atom)) if slot_atom[i] == A]
        assert len(own) == 4, (name, A, len(own))
        for j, i in enumerate(own):                                           # rep 1: plain OAOs on A's slots
            slots1_assign[i] = (["s", "p1", "p2", "p3"][j], None, cols[A][:, j])

    for A in sorted(heavy, key=lambda A: (len(heavy_nb[A]), -Rc[A, 0])):
        if A not in seen:
            dfs(A, None)
    for Hh in [A for A in range(len(sym)) if sym[A] == "H" and A not in h_owner]:
        emit_h(Hh)
    assert len(slot_atom) == n, (name, len(slot_atom), n)
    B0 = np.stack([v for _, _, v in slots0], 1)
    rep1 = []
    for i in range(n):
        if slots0[i][0] == "Hs":
            rep1.append(slots0[i])
        else:
            rep1.append(slots1_assign[i])
    B1 = np.stack([v for _, _, v in rep1], 1)
    Bs = []
    for Bk in (B0, B1):
        err = np.abs(Bk.T @ Bk - np.eye(n)).max()
        assert err < 1e-3, (name, err)
        u_, _, vt = np.linalg.svd(Bk)                                         # exact orthogonality (polar factor)
        Bs.append(u_ @ vt)
    B0, B1 = Bs
    t0 = {"Hs": 0, "parent": 1, "H": 2, "lp": 3, "child": 4}
    t1 = {"Hs": 0, "s": 6, "p1": 7, "p2": 8, "p3": 9}
    st = np.array([[t0[k] for k, _, _ in slots0], [t1[k] for k, _, _ in rep1]], np.int8)
    partner = np.array([pt if k in ("parent", "H", "child") else -1 for k, pt, _ in slots0], np.int16)
    sdir = np.zeros((n, 3), np.float32)                                      # rep-0 hybrid p-direction (canonical)
    for i, (k, _, _) in enumerate(slots0):
        if k != "Hs":
            A = slot_atom[i]
            coef = np.linalg.lstsq(cols[A], B0[:, i], rcond=None)[0]           # (s, p1, p2, p3) coefficients
            pn = np.linalg.norm(coef[1:])
            sdir[i] = coef[1:] / pn if pn > 1e-6 else 0.0
    return {"B": np.stack([B0, B1]), "slot_atom": np.array(slot_atom, np.int16), "slot_type": st,
            "slot_partner": partner, "slot_dir": sdir,
            "atom_Z": np.array([ZNUM[s] for s in sym], np.int8), "coords": Rc.astype(np.float32), "bond": bond}


_OPTS = {}


def _one(name):
    try:
        return name, build_frame(name, **_OPTS)
    except Exception as e:  # noqa: BLE001
        return name, repr(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-procs", type=int, default=24)
    ap.add_argument("--shapes", nargs="*", default=None, help="e.g. 16_13; default all")
    ap.add_argument("--out-dir", default="rhf_frames")
    ap.add_argument("--bond-mode", default="distance", choices=["distance", "mayer", "either"])
    ap.add_argument("--bo-thresh", type=float, default=0.3)
    args = ap.parse_args()
    _OPTS.update(bond_mode=args.bond_mode, bo_thresh=args.bo_thresh)
    idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
    groups = defaultdict(list)
    for k, (n, no, nv) in idx.items():
        groups[(no, nv)].append(k)
    out = ROOT / args.out_dir
    out.mkdir(exist_ok=True)
    t0 = time.time()
    with Pool(args.n_procs) as pool:
        for (no, nv), names in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            if args.shapes and f"{no}_{nv}" not in args.shapes:
                continue
            res = pool.map(_one, sorted(names), chunksize=4)
            ok = [(k, r) for k, r in res if isinstance(r, dict)]
            bad = [(k, r) for k, r in res if not isinstance(r, dict)]
            if not ok:
                print(f"  ({no},{nv}): all {len(bad)} failed, e.g. {bad[0]}", flush=True)
                continue
            na = max(len(r["atom_Z"]) for _, r in ok)
            M = len(ok)
            Z = np.zeros((M, na), np.int8)
            Cc = np.zeros((M, na, 3), np.float32)
            Bd = np.zeros((M, na, na), bool)
            for i, (_, r) in enumerate(ok):
                a = len(r["atom_Z"])
                Z[i, :a], Cc[i, :a], Bd[i, :a, :a] = r["atom_Z"], r["coords"], r["bond"]
            np.savez(out / f"{no}_{nv}.npz", names=np.array([k for k, _ in ok]),
                     B=np.stack([r["B"] for _, r in ok]).astype(np.float32),
                     slot_atom=np.stack([r["slot_atom"] for _, r in ok]),
                     slot_type=np.stack([r["slot_type"] for _, r in ok]),
                     slot_partner=np.stack([r["slot_partner"] for _, r in ok]),
                     slot_dir=np.stack([r["slot_dir"] for _, r in ok]),
                     atom_Z=Z, coords=Cc, bond=Bd, natm=np.array([len(r["atom_Z"]) for _, r in ok]))
            print(f"  ({no},{nv}) n={no+nv}: {M} ok, {len(bad)} failed{' e.g. ' + str(bad[0]) if bad else ''}"
                  f"  ({time.time()-t0:.0f}s)", flush=True)
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
