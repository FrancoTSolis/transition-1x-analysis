#!/usr/bin/env python3
"""Debug NaN in training."""
import torch
from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig
from pretrain.train import LUCJLoss

ds = CCSDAmplitudeDataset("jobs")

samples = []
seen_norbs = set()
for i in range(len(ds)):
    s = ds[i]
    if s["norb"] not in seen_norbs and len(samples) < 4:
        samples.append(s)
        seen_norbs.add(s["norb"])
        print(f"Sample {i}: norb={s['norb']}, nocc={s['nocc']}, nvirt={s['nvirt']}")
    if len(samples) >= 4:
        break

batch = CCSDAmplitudeDataset.collate_fn(samples)
print(f"\nBatch t2 shape: {batch['t2'].shape}")
print(f"norbs: {batch['norbs']}")

config = ModelConfig(embed_dim=128, num_layers=1, num_heads=4, n_reps=2)
model = PretrainingModel(config)

out = model(batch)
print(f"j_pred nan: {torch.isnan(out['j_pred']).any().item()}")
print(f"u_real_pred nan: {torch.isnan(out['u_real_pred']).any().item()}")

criterion = LUCJLoss()
loss, j_loss, u_loss = criterion(
    j_pred=out["j_pred"],
    j_target=batch["j_target"],
    u_real_pred=out["u_real_pred"],
    u_real_target=batch["u_real_target"],
    u_imag_pred=out["u_imag_pred"],
    u_imag_target=batch["u_imag_target"],
    norbs=batch["norbs"],
    max_norb=out["j_pred"].shape[-1],
)
print(f"loss={loss.item():.4f}, j={j_loss.item():.4f}, u={u_loss.item():.4f}")

# Now do a gradient step
loss.backward()
print("Backward OK")

# Check for nan in grads
for name, p in model.named_parameters():
    if p.grad is not None and torch.isnan(p.grad).any():
        print(f"NaN grad in {name}")
print("Grad check done")
