#!/usr/bin/env python3
"""Diagnose the U and J loss - inspect targets, predictions, and scales."""
import torch
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig

ds = CCSDAmplitudeDataset("jobs")

# Load a few samples of different sizes
indices = [0, 100, 500, 1000]
samples = [ds[i] for i in indices]
batch = CCSDAmplitudeDataset.collate_fn(samples)

print("=== TARGET STATISTICS ===")
for key in ["j_target", "u_real_target", "u_imag_target"]:
    t = batch[key]
    print(f"{key}: shape={t.shape}, min={t.min():.4f}, max={t.max():.4f}, "
          f"mean={t.mean():.4f}, std={t.std():.4f}, "
          f"fro_norm_per_sample={torch.sqrt((t**2).sum(dim=tuple(range(1, t.ndim)))).tolist()}")

print()

# Check: what's the MSE if we predict all zeros?
norbs = batch["norbs"]
B = len(norbs)
j_target = batch["j_target"]
u_real_target = batch["u_real_target"]
u_imag_target = batch["u_imag_target"]

j_zero_loss = 0.0
u_zero_loss = 0.0
for i in range(B):
    n = norbs[i].item()
    j_zero_loss += (j_target[i, :, :, :n, :n] ** 2).sum().item()
    u_zero_loss += (u_real_target[i, :, :n, :n] ** 2).sum().item()
    u_zero_loss += (u_imag_target[i, :, :n, :n] ** 2).sum().item()

print(f"MSE loss if predicting ALL ZEROS:")
print(f"  J: {j_zero_loss / B:.4f}")
print(f"  U: {u_zero_loss / B:.4f}")
print(f"  Total: {(j_zero_loss + u_zero_loss) / B:.4f}")
print()

# Check: what does U_real actually look like? Is it near-identity?
print("=== SAMPLE U TARGETS ===")
for i in range(min(2, B)):
    n = norbs[i].item()
    print(f"\nSample {indices[i]}: norb={n}, name={batch['names'][i]}")
    for rep in range(2):
        ur = u_real_target[i, rep, :n, :n]
        ui = u_imag_target[i, rep, :n, :n]
        print(f"  Rep {rep}:")
        print(f"    U_real: diag_mean={torch.diag(ur).mean():.4f}, diag_std={torch.diag(ur).std():.4f}, "
              f"off_diag_mean={(ur - torch.diag(torch.diag(ur))).abs().mean():.4f}")
        print(f"    U_imag: diag_mean={torch.diag(ui).mean():.4f}, diag_std={torch.diag(ui).std():.4f}, "
              f"off_diag_mean={(ui - torch.diag(torch.diag(ui))).abs().mean():.4f}")
        print(f"    U_real fro_norm={torch.norm(ur, 'fro'):.4f}, U_imag fro_norm={torch.norm(ui, 'fro'):.4f}")
        
        # Check if U is close to identity
        eye = torch.eye(n)
        print(f"    ||U_real - I||_F = {torch.norm(ur - eye, 'fro'):.4f}")

print()

# Now check model predictions
config = ModelConfig(embed_dim=128, num_layers=4, num_heads=8, n_reps=2, dropout=0.0)
model = PretrainingModel(config)
model.eval()

# Try loading best checkpoint
import glob
for ckpt_dir in ["checkpoints", "checkpoints_hf", "checkpoints_mo"]:
    ckpt_path = f"{ckpt_dir}/best.pt"
    try:
        ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        model.load_state_dict(ckpt["model_state_dict"])
        print(f"Loaded checkpoint from {ckpt_path}")
        break
    except:
        continue

with torch.no_grad():
    out = model(batch)

print("\n=== PREDICTION STATISTICS ===")
for key in ["j_pred", "u_real_pred", "u_imag_pred"]:
    t = out[key]
    print(f"{key}: min={t.min():.4f}, max={t.max():.4f}, mean={t.mean():.6f}, std={t.std():.4f}")

# Heatmaps
fig, axes = plt.subplots(4, 6, figsize=(30, 20))
i = 0  # first sample
n = norbs[i].item()
name = batch["names"][i]

# Row 0: J targets
for spin, label in enumerate(["J_aa", "J_ab"]):
    ax = axes[0, spin]
    im = ax.imshow(j_target[i, 0, spin, :n, :n].numpy(), cmap="RdBu_r", aspect="equal")
    ax.set_title(f"Target {label} (rep 0)")
    plt.colorbar(im, ax=ax, shrink=0.8)

# Row 0: J predictions
for spin, label in enumerate(["J_aa", "J_ab"]):
    ax = axes[0, spin + 2]
    im = ax.imshow(out["j_pred"][i, 0, spin, :n, :n].numpy(), cmap="RdBu_r", aspect="equal")
    ax.set_title(f"Pred {label} (rep 0)")
    plt.colorbar(im, ax=ax, shrink=0.8)

# Row 0: J diff
for spin, label in enumerate(["J_aa", "J_ab"]):
    ax = axes[0, spin + 4]
    diff = (out["j_pred"][i, 0, spin, :n, :n] - j_target[i, 0, spin, :n, :n]).numpy()
    im = ax.imshow(diff, cmap="RdBu_r", aspect="equal")
    ax.set_title(f"Diff {label}")
    plt.colorbar(im, ax=ax, shrink=0.8)

# Row 1: U_real target, pred, diff
for col, (data, title) in enumerate([
    (u_real_target[i, 0, :n, :n], "Target U_real"),
    (out["u_real_pred"][i, 0, :n, :n], "Pred U_real"),
    (out["u_real_pred"][i, 0, :n, :n] - u_real_target[i, 0, :n, :n], "Diff U_real"),
]):
    ax = axes[1, col]
    im = ax.imshow(data.numpy(), cmap="RdBu_r", aspect="equal")
    ax.set_title(f"{title} (rep 0)")
    plt.colorbar(im, ax=ax, shrink=0.8)

# Row 1: U_imag target, pred, diff
for col, (data, title) in enumerate([
    (u_imag_target[i, 0, :n, :n], "Target U_imag"),
    (out["u_imag_pred"][i, 0, :n, :n], "Pred U_imag"),
    (out["u_imag_pred"][i, 0, :n, :n] - u_imag_target[i, 0, :n, :n], "Diff U_imag"),
]):
    ax = axes[1, col + 3]
    im = ax.imshow(data.numpy(), cmap="RdBu_r", aspect="equal")
    ax.set_title(f"{title} (rep 0)")
    plt.colorbar(im, ax=ax, shrink=0.8)

# Row 2: Histograms of target values
ax = axes[2, 0]
ax.hist(j_target[i, 0, 0, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="J_aa")
ax.hist(j_target[i, 0, 1, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="J_ab")
ax.set_title("J target value distribution")
ax.legend()

ax = axes[2, 1]
ax.hist(u_real_target[i, 0, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="U_real")
ax.hist(u_imag_target[i, 0, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="U_imag")
ax.set_title("U target value distribution")
ax.legend()

ax = axes[2, 2]
ax.hist(out["j_pred"][i, 0, 0, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="J_aa pred")
ax.set_title("J pred value distribution")
ax.legend()

ax = axes[2, 3]
ax.hist(out["u_real_pred"][i, 0, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="U_real pred")
ax.hist(out["u_imag_pred"][i, 0, :n, :n].numpy().ravel(), bins=50, alpha=0.7, label="U_imag pred")
ax.set_title("U pred value distribution")
ax.legend()

# Row 3: diagonal analysis of U
ax = axes[2, 4]
diag_ur = torch.diag(u_real_target[i, 0, :n, :n]).numpy()
ax.plot(diag_ur, 'o-', markersize=3)
ax.axhline(0, color='gray', linestyle='--')
ax.set_title("U_real diagonal (target)")

ax = axes[2, 5]
diag_ui = torch.diag(u_imag_target[i, 0, :n, :n]).numpy()
ax.plot(diag_ui, 'o-', markersize=3)
ax.axhline(0, color='gray', linestyle='--')
ax.set_title("U_imag diagonal (target)")

# Hide unused subplots
for r in [3]:
    for c in range(6):
        axes[r, c].axis('off')

plt.suptitle(f"Diagnosis: {name} (norb={n})", fontsize=14)
plt.tight_layout()
plt.savefig("diagnose_loss.png", dpi=150)
print(f"\nSaved diagnose_loss.png")
