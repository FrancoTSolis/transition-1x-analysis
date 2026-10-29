"""Validate LUCJGradGPU against ffsim (complex128): energy, sigma = H psi, and the full gradient by central
finite differences of the ffsim energy.  Test Hamiltonians: leading norb x norb block of a TS molecule's
active space (a valid real Hamiltonian) with k electrons per spin.

    CUDA_VISIBLE_DEVICES=2 python -m ucj_init.test_lucj_grad --norb 8 --k 4
"""
from __future__ import annotations

import argparse
import time

import numpy as np

from ucj_init.lucj_grad_gpu import LUCJGradGPU, n_params, square_pairs, uz_from_x, x_from_uz


def ffsim_energy_x(x, h1, eri, const, n, k, t1):
    import ffsim
    from pretrain.rl.energy import make_ucj_op
    U, Z = uz_from_x(x, n, 2)
    op = make_ucj_op(Z, U, "square", t1=t1)
    ham = ffsim.MolecularHamiltonian(h1, eri, const)
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(n, (k, k)), op, norb=n, nelec=(k, k))
    return float(np.vdot(psi, ffsim.linear_operator(ham, norb=n, nelec=(k, k)) @ psi).real), psi, ham


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="C2H3N_rxn2858_TS")
    ap.add_argument("--norb", type=int, default=8)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--scale", type=float, default=0.3)
    ap.add_argument("--fd", type=int, default=40, help="number of gradient components checked by FD")
    args = ap.parse_args()
    import ffsim
    from ucj_init.common import ROOT
    d = np.load(ROOT / "rhf_hamiltonians" / f"{args.name}.npz")
    n, k = args.norb, args.k
    h1 = d["one_body"][:n, :n].astype(np.float64)
    eri = d["two_body"][:n, :n, :n, :n].astype(np.float64)
    const = float(d["constant"])
    rng = np.random.default_rng(1)
    t1 = 0.05 * rng.normal(size=(k, n - k))
    x = args.scale * rng.normal(size=n_params(n, 2))
    # round-trip of the parameterization against ffsim's
    from pretrain.rl.energy import make_ucj_op
    U, Z = uz_from_x(x, n, 2)
    op = make_ucj_op(Z, U, "square")
    xf = op.to_parameters(interaction_pairs=square_pairs(n))
    print("x vs ffsim.to_parameters (principal log; equal when |K| small):", np.abs(xf - x).max(),
          " x_from_uz:", np.abs(x_from_uz(U, Z) - x).max())

    t0 = time.time()
    eng = LUCJGradGPU(h1, eri, const, n, (k, k), t1=t1, verbose=True)
    print(f"build {time.time() - t0:.2f}s  Ha asym {eng.Ha_asym:.1e}")
    E, g = eng.energy_grad(x)
    Ef, psi, ham = ffsim_energy_x(x, h1, eri, const, n, k, t1)
    print(f"E gpu {E:.10f}  ffsim {Ef:.10f}  diff {E - Ef:.2e}")
    # sigma check on the ffsim state
    import torch
    P = torch.as_tensor(psi.reshape(eng.dim, eng.dim), device=eng.device, dtype=torch.complex64)
    S = torch.empty_like(P)
    eng.sigma(P, S)
    sig_ref = ffsim.linear_operator(ham, norb=n, nelec=(k, k)) @ psi - (const + eng.E0) * psi
    err = np.abs(S.cpu().numpy().reshape(-1) - sig_ref).max()
    print(f"sigma max abs err {err:.2e}  (max |sigma| {np.abs(sig_ref).max():.2e})")
    # gradient by central FD of ffsim energies on a subset of components
    idx = rng.choice(len(x), size=min(args.fd, len(x)), replace=False)
    h = 1e-5
    fd = []
    for i in idx:
        xp = x.copy(); xp[i] += h
        xm = x.copy(); xm[i] -= h
        fd.append((ffsim_energy_x(xp, h1, eri, const, n, k, t1)[0] - ffsim_energy_x(xm, h1, eri, const, n, k, t1)[0]) / (2 * h))
    fd = np.array(fd)
    print(f"grad: max |g - fd| {np.abs(g[idx] - fd).max():.2e}   max |fd| {np.abs(fd).max():.2e}   "
          f"rel {np.linalg.norm(g[idx] - fd) / np.linalg.norm(fd):.2e}")
    for j, i in enumerate(idx[:8]):
        print(f"   {i:4d}  gpu {g[i]:+.8f}  fd {fd[j]:+.8f}")
    t0 = time.time()
    for _ in range(3):
        eng.energy_grad(x)
    print(f"energy+grad {(time.time() - t0) / 3:.3f}s")


if __name__ == "__main__":
    main()
