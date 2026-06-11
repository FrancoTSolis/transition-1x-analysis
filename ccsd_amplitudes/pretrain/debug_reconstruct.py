#!/usr/bin/env python3
"""Debug: why is LUCJ reconstruction residual ~0.96?

Compare reconstruction quality across DF settings to (a) validate the
einsum formula and (b) see how reconstruction depends on n_reps / sparsity.
"""
from pathlib import Path
import numpy as np


def recon(Z, U, nocc):
    return 1j * np.einsum("kpq,kap,kip,kbq,kjq->ijab",
                          Z, U, U.conj(), U, U.conj())[:nocc, :nocc, nocc:, nocc:]


def rel(t2_rec, t2):
    return np.linalg.norm(t2_rec.real - t2) / np.linalg.norm(t2)


def main():
    import ffsim
    s = np.load(sorted(Path("rhf_cache").glob("*.npz"))[0])
    t1, t2 = s["t1"], s["t2"]
    nocc, nvirt, norb = int(s["nocc"]), int(s["nvirt"]), int(s["norb"])
    print(f"norb={norb}, nocc={nocc}, nvirt={nvirt}, ||t2||={np.linalg.norm(t2):.3f}\n")

    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    pairs_ab = [(p, p) for p in range(norb)]

    configs = [
        ("exact DF (no trunc, no constraints)", dict(n_reps=None, optimize=False, interaction_pairs=None)),
        ("n_reps=2, no constraints, no opt",     dict(n_reps=2, optimize=False, interaction_pairs=None)),
        ("n_reps=2, no constraints, opt",        dict(n_reps=2, optimize=True, interaction_pairs=None)),
        ("n_reps=2, LUCJ constraints, opt",      dict(n_reps=2, optimize=True, interaction_pairs=(pairs_aa, pairs_ab))),
        ("n_reps=8, no constraints, opt",        dict(n_reps=8, optimize=True, interaction_pairs=None)),
    ]
    for label, kw in configs:
        opt = kw.pop("optimize")
        try:
            extra = dict(optimize=opt)
            if opt:
                extra["options"] = dict(maxiter=100)
            op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2=t2, **kw, **extra)
            U = op.orbital_rotations
            Z = op.diag_coulomb_mats[:, 0]
            r = rel(recon(Z, U, nocc), t2)
            print(f"  {label:42s} n_reps_actual={U.shape[0]:2d}  resid={r:.4f}")
        except Exception as e:
            print(f"  {label:42s} ERROR: {e}")


if __name__ == "__main__":
    main()
