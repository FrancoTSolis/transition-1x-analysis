#!/usr/bin/env python3
"""Sanity: does the LUCJ (U,Z) reconstruct t2, and is U real?

Establishes the reconstruction-error FLOOR (best a model could achieve at
n_reps=2) and the zero-prediction baseline (the ceiling). If the floor is
well below the ceiling, 'reconstruct t2' is a meaningful, well-posed target.
"""
import sys
from pathlib import Path
import numpy as np


def main():
    import ffsim
    from scipy.linalg import expm

    def reconstruct_t2(Z, U, nocc):
        return 1j * np.einsum("kpq,kap,kip,kbq,kjq->ijab",
                              Z, U, U.conj(), U, U.conj())[:nocc, :nocc, nocc:, nocc:]

    files = sorted(Path("rhf_cache").glob("*.npz"))[:8]
    print(f"Checking {len(files)} molecules\n")
    floors, ceils, imag_frac = [], [], []
    for f in files:
        s = np.load(f)
        t1, t2 = s["t1"], s["t2"]
        nocc, nvirt, norb = int(s["nocc"]), int(s["nvirt"]), int(s["norb"])

        pairs_aa = [(p, p + 1) for p in range(norb - 1)]
        pairs_ab = [(p, p) for p in range(norb)]
        op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
            t2=t2, t1=t1, n_reps=2, interaction_pairs=(pairs_aa, pairs_ab),
            optimize=True, options=dict(maxiter=50),
        )
        U = op.orbital_rotations            # (2, norb, norb) complex
        Z = op.diag_coulomb_mats[:, 0]      # (2, norb, norb) real
        t2_rec = reconstruct_t2(Z, U, nocc)

        resid = np.linalg.norm(t2_rec.real - t2) / np.linalg.norm(t2)
        zero = np.linalg.norm(t2) / np.linalg.norm(t2)  # =1 (predict zeros)
        floors.append(resid)
        ceils.append(zero)
        imag_frac.append(np.linalg.norm(U.imag) / np.linalg.norm(U.real))

        print(f"  {f.stem[:28]:28s} recon_resid={resid:.3f}  U_imag/U_real={imag_frac[-1]:.3f}")

    print(f"\nFLOOR  (LUCJ n_reps=2 reconstruction error): mean={np.mean(floors):.3f}")
    print(f"CEILING(predict zeros, relative error)     : mean={np.mean(ceils):.3f}")
    print(f"U imag/real magnitude ratio                : mean={np.mean(imag_frac):.3f}")
    print(f"\n=> reconstruction target spans [{np.mean(floors):.2f}, 1.00]; "
          f"a model beating ~1.0 is learning real structure.")
    if np.mean(imag_frac) < 0.05:
        print("=> U is effectively REAL orthogonal: can parametrize via expm(antisym).")
    else:
        print("=> U has non-trivial imaginary part: parametrize via expm(anti-Hermitian).")


if __name__ == "__main__":
    main()
