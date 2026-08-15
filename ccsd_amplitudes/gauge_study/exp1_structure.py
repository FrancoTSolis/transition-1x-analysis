#!/usr/bin/env python3
"""Experiment 1: verify the mathematical structure of the stored labels.

Checks, on real rhf_dataset_n29 molecules + their rhf_targets:
  A. numpy reimplementation of the exact DF reproduces stored (Z, U):
     Z should match exactly; U should match ONLY up to gauge
     (per-column phase + permutation) -> quantifies how much of the stored
     label is pure gauge picked by LAPACK.
  B. structural identities: Z_1 = -Z_0, rank(Z_0) = 1, U_1 ~ conj(U_0).
  C. the invariant object: reconstruct_t2(Z, U) equals the best rank-1
     approximation of the t2 matrix (lam0 * v0 v0^T).

Run from ccsd_amplitudes/:  python3 -m gauge_study.exp1_structure
"""
import numpy as np

from .common import (align_unitaries, exact_df_t2, invariant_B, load_n29,
                     reconstruct_t2)


def main():
    rows = load_n29(max_mols=50)
    print(f"{len(rows)} molecules (norb=29, nocc=17, nvirt=12)\n")

    dZ_exact, dU_raw, dU_aligned = [], [], []
    dZ_antisym, Z_rank1_resid, dU_conj = [], [], []
    dB = []
    for r in rows:
        ref = exact_df_t2(r["t2"], n_reps=2)

        # A. reproduce stored targets
        dZ_exact.append(np.abs(ref["Z"] - r["Z"]).max())
        for k in range(2):
            _, _, d_al = align_unitaries(r["U"][k], ref["U"][k])
            dU_aligned.append(d_al / np.linalg.norm(r["U"][k]))
            dU_raw.append(
                np.linalg.norm(r["U"][k] - ref["U"][k]) / np.linalg.norm(r["U"][k]))

        # B. structure
        dZ_antisym.append(np.abs(r["Z"][0] + r["Z"][1]).max())
        s = np.linalg.svd(r["Z"][0], compute_uv=False)
        Z_rank1_resid.append(s[1] / s[0] if s[0] > 0 else 0.0)
        _, _, d_conj = align_unitaries(r["U"][1], r["U"][0].conj())
        dU_conj.append(d_conj / np.linalg.norm(r["U"][1]))

        # C. invariant object
        t2_hat = reconstruct_t2(r["Z"], r["U"], r["nocc"])
        B = invariant_B(r["t2"])
        dB.append(np.linalg.norm(t2_hat - B) / np.linalg.norm(B))

    def stat(name, arr):
        arr = np.array(arr)
        print(f"  {name:44s} median={np.median(arr):.2e}  max={arr.max():.2e}")

    print("A. reimplementation vs stored targets")
    stat("max|Z_reimpl - Z_stored|", dZ_exact)
    stat("||U_reimpl - U_stored|| / ||U||  (RAW)", dU_raw)
    stat("||U_reimpl - U_stored|| / ||U||  (ALIGNED)", dU_aligned)

    print("\nB. structural identities of the stored labels")
    stat("max|Z_1 + Z_0|            (Z_1 = -Z_0 ?)", dZ_antisym)
    stat("sigma_2/sigma_1 of Z_0    (rank-1 ?)", Z_rank1_resid)
    stat("aligned ||U_1 - conj(U_0)|| / ||U||", dU_conj)

    print("\nC. invariant content = best rank-1 approx of t2 matrix")
    stat("||t2hat(Z,U) - lam0 v0 v0^T|| / ||B||", dB)


if __name__ == "__main__":
    main()
