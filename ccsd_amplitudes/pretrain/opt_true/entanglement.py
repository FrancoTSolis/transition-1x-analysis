#!/usr/bin/env python3
"""How MPS-friendly is the LUCJ state in the chemistry frame vs the MO basis?

For a small molecule (exact statevector), build the LUCJ state psi from given (U, Z) parameters, express it in an
orbital basis (MO, or the chemistry frame B of rep 0 / rep 1, chain-ordered), and compute the Schmidt spectrum
across every orbital cut of a 1D ordering:
  interleaved  ... a_p b_p a_{p+1} b_{p+1} ...   (spin orbitals of one spatial orbital adjacent)
  blocked      a_0 .. a_{n-1} b_0 .. b_{n-1}     (qiskit / ffsim JW qubit order used by CircuitMPS)
Reported per cut: von Neumann entropy and the bond dimension chi needed to keep 1 - 1e-3 / 1 - 1e-4 of the norm.
The fermionic reordering signs between orderings are constant within each particle-number block, so the
singular values are exact.

Usage: python3 -m pretrain.opt_true.entanglement --name C2H3N_rxn2857_P --params label
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def occ_strings(n, k):
    import math
    import ffsim
    addrs = np.arange(math.comb(n, k))
    out = ffsim.addresses_to_strings(addrs, n, k, bitstring_type=ffsim.BitstringType.INT)
    return np.asarray(out, dtype=np.int64)


def _sv(M):
    """Singular values via the eigenvalues of the smaller Gram matrix (much faster than a full SVD)."""
    G = M @ M.conj().T if M.shape[0] <= M.shape[1] else M.conj().T @ M
    w = np.linalg.eigvalsh(G)
    return np.sqrt(np.clip(w, 0, None))


def schmidt_interleaved(psi, n, na, nb, cut):
    """Singular values across the cut: left = spatial orbitals [0, cut) (both spins)."""
    sa, sb = occ_strings(n, na), occ_strings(n, nb)
    lm = (1 << cut) - 1
    aL, aR = sa & lm, sa >> cut
    bL, bR = sb & lm, sb >> cut
    pc = np.vectorize(lambda x: bin(int(x)).count("1"))
    naL, nbL = pc(aL), pc(bL)
    svals = []
    for x in range(min(na, cut) + 1):
        ia = np.where(naL == x)[0]
        if not len(ia):
            continue
        for y in range(min(nb, cut) + 1):
            ib = np.where(nbL == y)[0]
            if not len(ib):
                continue
            # block rows: (aL, bL) pairs, cols: (aR, bR) pairs
            uaL, ra = np.unique(aL[ia], return_inverse=True)
            uaR, ca = np.unique(aR[ia], return_inverse=True)
            ubL, rb = np.unique(bL[ib], return_inverse=True)
            ubR, cb = np.unique(bR[ib], return_inverse=True)
            M = np.zeros((len(uaL) * len(ubL), len(uaR) * len(ubR)), dtype=psi.dtype)
            sub = psi[np.ix_(ia, ib)]
            rows = (ra[:, None] * len(ubL) + rb[None, :])
            cols = (ca[:, None] * len(ubR) + cb[None, :])
            M[rows, cols] = sub
            svals.append(_sv(M))
    return np.sort(np.concatenate(svals))[::-1]


def schmidt_blocked(psi, n, na, nb, cut):
    """Blocked order a_0..a_{n-1} b_0..b_{n-1}; cut counted in spin orbitals (0..2n)."""
    if cut <= n:
        sa = occ_strings(n, na)
        lm = (1 << cut) - 1
        aL, aR = sa & lm, sa >> cut
        pc = np.vectorize(lambda x: bin(int(x)).count("1"))
        naL = pc(aL)
        svals = []
        for x in range(min(na, cut) + 1):
            ia = np.where(naL == x)[0]
            if not len(ia):
                continue
            uaL, ra = np.unique(aL[ia], return_inverse=True)
            uaR, ca = np.unique(aR[ia], return_inverse=True)
            # rows aL, cols (aR, b)
            M = np.zeros((len(uaL), len(uaR) * psi.shape[1]), dtype=psi.dtype)
            for r, c, i in zip(ra, ca, ia):
                M[r, c * psi.shape[1]:(c + 1) * psi.shape[1]] = psi[i]
            svals.append(_sv(M))
        return np.sort(np.concatenate(svals))[::-1]
    return schmidt_blocked(psi.T.copy(), n, nb, na, 2 * n - cut)


def stats(s):
    p = s ** 2
    p = p / p.sum()
    ent = float(-(p[p > 1e-16] * np.log(p[p > 1e-16])).sum())
    c = np.cumsum(p)
    return {"S": round(ent, 3), "chi_1e-3": int(np.searchsorted(c, 1 - 1e-3) + 1), "chi_1e-4": int(np.searchsorted(c, 1 - 1e-4) + 1)}


def main():
    import ffsim
    import torch
    from pretrain.rl.energy import make_ucj_op
    from pretrain.rl.hamiltonian import load_hamiltonian
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="C2H3N_rxn2857_P")
    ap.add_argument("--params", default="label", choices=["label", "frame_adam"])
    ap.add_argument("--labels-dir", default="rhf_targets_compressed_small")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    name = args.name
    ham, n, nelec, e_hf, e_ccsd = load_hamiltonian(ROOT / "rhf_hamiltonians", name)
    na, nb = nelec
    d = np.load(ROOT / "rhf_dataset" / f"{name}.npz")
    t1 = d["t1"].astype(np.float64)
    if args.params == "label":
        lab = np.load(ROOT / args.labels_dir / "square_reg0.005" / f"{name}.npz")
        U, Z = lab["U_re"] + 1j * lab["U_im"], lab["Z"]
    else:
        from pretrain.opt_true.eval_slot_energy import in_frame_adam
        from pretrain.opt_true.train_slot_all import Bucket
        idx = json.load(open(ROOT / "rhf_dataset" / "_index.json"))
        _, no, nv = idx[name]
        x = Bucket(no, nv, [name]).batch([0], "cpu")
        from pretrain.opt_true import dftorch as D
        Uo, Zo = in_frame_adam(x["U0"].to(torch.complex128), x["t2"].double(), D.square_mask(n, "cpu", torch.float64),
                               0.005, x["zref"].double(), 300, nocc=no)
        U, Z = Uo[0].numpy(), Zo[0].numpy()
    W, _, Vh = np.linalg.svd(U)
    U = W @ Vh
    op = make_ucj_op(Z, U, "square", t1=t1)
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(n, nelec), op, norb=n, nelec=nelec)
    E = float(np.vdot(psi, ffsim.linear_operator(ham, norb=n, nelec=nelec) @ psi).real)
    print(f"{name} ({args.params}): corr% {(e_hf - E) / (e_hf - e_ccsd) * 100:.1f}", flush=True)
    import math
    psi = psi.reshape(math.comb(n, na), math.comb(n, nb))
    f = np.load(ROOT / "rhf_frames" / f"{na}_{n - na}.npz")
    i = [str(k) for k in f["names"]].index(name)
    bases = {"MO (energy order)": np.eye(n), "frame rep0 (hybrids)": f["B"][i, 0].astype(np.float64),
             "frame rep1 (OAOs)": f["B"][i, 1].astype(np.float64)}
    res = {}
    for tag, B in bases.items():
        u_, _, vt = np.linalg.svd(B)
        B = u_ @ vt
        v = ffsim.apply_orbital_rotation(psi.reshape(-1), B.conj().T, norb=n, nelec=nelec).reshape(psi.shape)
        inter = [stats(schmidt_interleaved(v, n, na, nb, c)) for c in range(1, n)]
        mid = n // 2
        blk = stats(schmidt_blocked(v, n, na, nb, n))          # the alpha | beta cut of the blocked order
        res[tag] = {"interleaved": inter, "blocked_mid": blk}
        mx = max(inter, key=lambda s: s["chi_1e-4"])
        print(f"  {tag:22s} interleaved: max S {max(s['S'] for s in inter):.2f}, max chi(1e-3) "
              f"{max(s['chi_1e-3'] for s in inter)}, max chi(1e-4) {mx['chi_1e-4']}, mid-cut {inter[mid - 1]} | "
              f"blocked alpha|beta cut: {blk}", flush=True)
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
