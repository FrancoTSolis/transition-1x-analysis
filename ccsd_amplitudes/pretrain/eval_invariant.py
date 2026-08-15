#!/usr/bin/env python3
"""Evaluate the deck_v3 invariant-retargeted model, vs the old kappa baseline.

Both models are scored on the SAME gauge-aware quantity: the invariant content
of the predicted LUCJ factorization,

    t2hat = i sum_k,pq Z_k[pq] U_k[ap] U_k*[ip] U_k[bq] U_k*[jq]
    (for exact n_reps=2 labels this equals B = lam0 v0 v0^T)

  new model:  (lam_hat, v_hat) -> quadratures -> eigh -> (U, Z) -> t2hat
  baseline :  (Z_pred, kappa_pred) -> U = expm(kappa) -> t2hat

reported as ||t2hat_pred - B_label|| / ||B_label||, plus direction cosine,
lam0 R^2, Z MSE, and gap-stratified medians.

Usage (from ccsd_amplitudes/):
    python3 -m pretrain.eval_invariant \
        --checkpoint checkpoints_invariant/best.pt \
        --baseline checkpoints_pretrain/best.pt [--max-mols 1000]
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset

from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import ModelConfig, PretrainingModel


def quadrature_np(M: np.ndarray, sign: int) -> np.ndarray:
    return 0.5 * (1 - sign * 1j) * (M + sign * 1j * M.T)


def uz_from_lam_v(lam: float, v: np.ndarray, nocc: int, nvirt: int):
    """ffsim's exact-DF construction from the invariant statistic."""
    norb = nocc + nvirt
    v = v / max(np.linalg.norm(v), 1e-12)
    M = np.zeros((norb, norb))
    M[nocc:, :nocc] = v.reshape(nocc, nvirt).T
    Zs, Us = [], []
    for sign, coeff in ((1, lam), (-1, -lam)):
        w, U = np.linalg.eigh(quadrature_np(M, sign))
        Zs.append(coeff * np.outer(w, w))
        Us.append(U)
    return np.array(Zs), np.array(Us)


def t2hat_from_uz(Z: np.ndarray, U: np.ndarray, nocc: int) -> np.ndarray:
    full = 1j * np.einsum(
        "kpq,kap,kip,kbq,kjq->ijab", Z, U, U.conj(), U, U.conj(), optimize=True)
    return full[:nocc, :nocc, nocc:, nocc:].real


def load_model(path: str, invariant: bool, device, args) -> PretrainingModel:
    cfg = ModelConfig(
        embed_dim=args.embed_dim, num_layers=args.num_layers,
        num_heads=args.num_heads, n_reps=2, dropout=0.0,
        predict_invariant=invariant,
    )
    model = PretrainingModel(cfg).to(device)
    sd = torch.load(path, map_location=device, weights_only=False)["model_state_dict"]
    missing, unexpected = model.load_state_dict(sd, strict=False)
    if missing or unexpected:
        print(f"  [{Path(path).name}] missing={len(missing)} unexpected={len(unexpected)}")
    model.eval()
    return model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True,
                    help="invariant-mode checkpoint (deck_v3)")
    ap.add_argument("--baseline", default=None,
                    help="old supervised (Z,kappa) checkpoint for comparison")
    ap.add_argument("--data-dir", default="rhf_dataset")
    ap.add_argument("--targets-dir", default="rhf_targets")
    ap.add_argument("--inv-targets-dir", default="rhf_inv_targets")
    ap.add_argument("--val-split", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--batch-size", type=int, default=12)
    ap.add_argument("--max-mols", type=int, default=0, help="0 = full val set")
    ap.add_argument("--embed-dim", type=int, default=192)
    ap.add_argument("--num-layers", type=int, default=6)
    ap.add_argument("--num-heads", type=int, default=8)
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dataset = CCSDAmplitudeDataset(
        args.data_dir, targets_dir=args.targets_dir,
        inv_targets_dir=args.inv_targets_dir)
    gen = torch.Generator().manual_seed(args.seed)
    indices = torch.randperm(len(dataset), generator=gen).tolist()
    n_val = int(len(dataset) * args.val_split)
    val_indices = indices[len(dataset) - n_val:]
    if args.max_mols:
        val_indices = val_indices[:args.max_mols]
    loader = DataLoader(
        Subset(dataset, val_indices), batch_size=args.batch_size, shuffle=False,
        collate_fn=dataset.collate_fn, num_workers=4)
    print(f"val set: {len(val_indices)} molecules on {device}")

    model_new = load_model(args.checkpoint, invariant=True, device=device, args=args)
    model_old = (load_model(args.baseline, invariant=False, device=device, args=args)
                 if args.baseline else None)

    rows = []
    t0 = time.time()
    with torch.no_grad():
        for batch in loader:
            batch_dev = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                         for k, v in batch.items()}
            out_new = model_new(batch_dev)
            out_old = model_old(batch_dev) if model_old is not None else None
            B = batch_dev["t2"].shape[0]
            max_nocc = batch_dev["max_nocc"]
            for i in range(B):
                no = int(batch_dev["noccs"][i]); nv = int(batch_dev["nvirts"][i])
                n = no + nv
                lam0 = float(batch_dev["lam_target"][i])
                v0 = batch_dev["v_target"][i][:no, :nv].cpu().numpy().astype(np.float64)
                gap = float(batch_dev["gap"][i])
                B_label = lam0 * np.outer(v0.ravel(), v0.ravel())
                nrmB = np.linalg.norm(B_label)

                # ---- new model ----
                m_hat = out_new["v_pred"][i][:no, max_nocc:max_nocc + nv]
                m_hat = m_hat.cpu().numpy().astype(np.float64)
                lam_h = float(out_new["lam_pred"][i])
                vh = m_hat.ravel() / max(np.linalg.norm(m_hat), 1e-12)
                cos = abs(float(vh @ v0.ravel()))
                B_new = lam_h * np.outer(vh, vh)
                err_new = np.linalg.norm(B_new - B_label) / nrmB

                # ---- baseline: (Z, kappa) -> U -> t2hat ----
                err_old = None
                if out_old is not None:
                    idx = np.r_[np.arange(no), max_nocc + np.arange(nv)]
                    Zp = out_old["j_pred"][i, :, 0].cpu().numpy().astype(np.float64)
                    Zp = Zp[:, idx][:, :, idx]
                    kr = out_old["kappa_real_pred"][i].cpu().numpy().astype(np.float64)
                    ki = out_old["kappa_imag_pred"][i].cpu().numpy().astype(np.float64)
                    kr = kr[:, idx][:, :, idx]; ki = ki[:, idx][:, :, idx]
                    from scipy.linalg import expm
                    U = np.stack([expm(kr[k] + 1j * ki[k]) for k in range(2)])
                    t2h = t2hat_from_uz(Zp, U, no)
                    Bmat_old = t2h.transpose(0, 2, 1, 3).reshape(no * nv, no * nv)
                    err_old = np.linalg.norm(Bmat_old - B_label) / nrmB

                rows.append(dict(cos=cos, err_new=err_new, err_old=err_old,
                                 lam_h=lam_h, lam0=lam0, gap=gap, norb=n))

    cos = np.array([r["cos"] for r in rows])
    err_new = np.array([r["err_new"] for r in rows])
    lam_h = np.array([r["lam_h"] for r in rows])
    lam0 = np.array([r["lam0"] for r in rows])
    gap = np.array([r["gap"] for r in rows])
    r2 = 1 - ((lam_h - lam0) ** 2).sum() / ((lam0 - lam0.mean()) ** 2).sum()

    print(f"\n=== invariant model ({len(rows)} val molecules, "
          f"{time.time()-t0:.0f}s) ===")
    print(f"B rel err : median={np.median(err_new):.3f}  mean={err_new.mean():.3f}")
    print(f"|cos(v0)| : median={np.median(cos):.3f}  mean={cos.mean():.3f}")
    print(f"lam0      : R^2={r2:.3f}")
    if rows[0]["err_old"] is not None:
        err_old = np.array([r["err_old"] for r in rows])
        print(f"\n=== baseline (Z,kappa)->expm->t2hat, same metric ===")
        print(f"B rel err : median={np.median(err_old):.3f}  mean={err_old.mean():.3f}")

    print("\n--- stratified by outer gap (invariant model) ---")
    for lo, hi in [(0.0, 0.05), (0.05, 0.2), (0.2, 1.01)]:
        m = (gap >= lo) & (gap < hi)
        if m.sum() < 5:
            continue
        print(f"  gap [{lo},{hi}):  n={m.sum():5d}  "
              f"B err median={np.median(err_new[m]):.3f}  "
              f"cos median={np.median(cos[m]):.3f}")


if __name__ == "__main__":
    main()
