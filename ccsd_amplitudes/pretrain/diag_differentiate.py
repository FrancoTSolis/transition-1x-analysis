#!/usr/bin/env python3
"""Diagnose whether the model produces molecule-DISTINCT (U,Z).

Loads a checkpoint, runs the model on several molecules, and reports how much
the predicted kappa/Z vary across molecules vs how much t2 varies. If the
outputs barely vary while t2 varies a lot, the molecule-specific pathway
(transformer -> heads) is not contributing (the shared base dominates).
"""
import sys
from pathlib import Path

import numpy as np
import torch

from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig


def main():
    ckpt = sys.argv[1] if len(sys.argv) > 1 else "/tmp/of_ckpt4/last.pt"
    data_dir = sys.argv[2] if len(sys.argv) > 2 else "rhf_overfit"
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    cfg = ModelConfig(embed_dim=192, num_layers=4, num_heads=6, n_reps=2)
    model = PretrainingModel(cfg).to(dev)
    state = torch.load(ckpt, map_location=dev, weights_only=False)
    model.load_state_dict(state["model_state_dict"])
    model.eval()

    ds = CCSDAmplitudeDataset(data_dir, n_reps=2)
    batch = ds.collate_fn([ds[i] for i in range(min(8, len(ds)))])
    batch = {k: v.to(dev) if isinstance(v, torch.Tensor) else v for k, v in batch.items()}

    with torch.no_grad():
        out = model(batch)
    kr = out["kappa_real_pred"]   # (B, reps, N, N)
    jz = out["j_pred"][:, :, 0]   # (B, reps, N, N)
    t2 = batch["t2"]              # (B, ...)

    B = kr.shape[0]
    # variation across molecules vs mean magnitude
    def variation(x):
        xf = x.reshape(B, -1)
        across = (xf - xf.mean(0, keepdim=True)).norm(dim=1).mean().item()
        mag = xf.norm(dim=1).mean().item()
        return across, mag

    a_kr, m_kr = variation(kr)
    a_jz, m_jz = variation(jz)
    a_t2, m_t2 = variation(t2)

    print(f"Checkpoint: {ckpt}   molecules: {B}\n")
    print(f"{'quantity':>10} {'cross-mol var':>14} {'mean magnitude':>15} {'var/mag':>9}")
    for name, (a, m) in [("t2", (a_t2, m_t2)), ("kappa_re", (a_kr, m_kr)), ("Z", (a_jz, m_jz))]:
        print(f"{name:>10} {a:>14.4f} {m:>15.4f} {a/max(m,1e-8):>9.3f}")

    print("\nPairwise output distances (should be >0 if model differentiates):")
    krf = kr.reshape(B, -1)
    for i in range(min(4, B)):
        for j in range(i + 1, min(4, B)):
            dk = (krf[i] - krf[j]).norm().item()
            dt = (t2.reshape(B, -1)[i] - t2.reshape(B, -1)[j]).norm().item()
            print(f"  mol {i}-{j}:  ||Δkappa||={dk:.4f}   ||Δt2||={dt:.4f}")


if __name__ == "__main__":
    main()
