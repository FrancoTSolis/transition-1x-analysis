#!/usr/bin/env python3
"""Evaluate the J/Z prediction quality of a trained LUCJ model.

Reports R^2 (vs predict-zero baseline) for the diagonal-Coulomb Z, overall and
on the structurally-nonzero entries, plus relative Frobenius error. Confirms the
J/Z head is genuinely learning (not just exploiting Z's small scale).
"""
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset

from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig


def main():
    ckpt = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sup_smoke/best.pt"
    data_dir = sys.argv[2] if len(sys.argv) > 2 else "rhf_dataset_n29"
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    ds = CCSDAmplitudeDataset(data_dir, n_reps=2, targets_dir="rhf_targets")
    n_val = int(len(ds) * 0.1)
    g = torch.Generator().manual_seed(42)
    idx = torch.randperm(len(ds), generator=g).tolist()
    val = Subset(ds, idx[len(ds) - n_val:])
    loader = DataLoader(val, batch_size=8, shuffle=False, collate_fn=ds.collate_fn, num_workers=0)

    cfg = ModelConfig(embed_dim=192, num_layers=6, num_heads=8, n_reps=2)
    model = PretrainingModel(cfg).to(dev)
    model.load_state_dict(torch.load(ckpt, map_location=dev, weights_only=False)["model_state_dict"])
    model.eval()

    se = 0.0; base = 0.0; n = 0           # overall
    se_nz = 0.0; base_nz = 0.0            # nonzero-target entries
    rel_errs = []
    with torch.no_grad():
        for batch in loader:
            nb = {k: (v.to(dev) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
            out = model(nb)
            jp = out["j_pred"][:, :, 0]   # (B, reps, N, N)
            for i in range(jp.shape[0]):
                no = int(nb["noccs"][i]); nv = int(nb["nvirts"][i]); m = no + nv
                gi = torch.cat([torch.arange(no, device=dev), nb["max_nocc"] + torch.arange(nv, device=dev)])
                Zp = jp[i][:, gi][:, :, gi].cpu().numpy()
                Zt = nb["z_target"][i][:, :m, :m].cpu().numpy()
                se += ((Zp - Zt) ** 2).sum(); base += (Zt ** 2).sum(); n += Zt.size
                mask = np.abs(Zt) > 1e-6
                se_nz += ((Zp - Zt) ** 2 * mask).sum(); base_nz += (Zt ** 2 * mask).sum()
                rel_errs.append(np.linalg.norm(Zp - Zt) / max(np.linalg.norm(Zt), 1e-8))

    print(f"checkpoint: {ckpt}   val molecules: {len(val)}")
    print(f"Z overall : R^2 = {1 - se/base:.4f}   (MSE {se/n:.3e} vs baseline {base/n:.3e})")
    print(f"Z nonzero : R^2 = {1 - se_nz/base_nz:.4f}   ({int(mask.sum())}-ish nz entries/sample)")
    print(f"Z rel.err : mean ||Zp-Z*||/||Z*|| = {np.mean(rel_errs):.4f}")


if __name__ == "__main__":
    main()
