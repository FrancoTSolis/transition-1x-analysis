#!/usr/bin/env python3
"""Result panels for the deck, computed from the trained invariant model.

Runs the model over the held-out validation molecules once, caches the
per-molecule numbers, and renders two deliberately plain panels:

  parity.png     predicted vs exact interaction strength (one dot per molecule)
  agreement.png  distribution of agreement with the exact setup

Run from ccsd_amplitudes/ with the training venv:
    pretrain/.train_venv/bin/python3 ../gates_slides/make_results.py
"""
import sys
from pathlib import Path

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ccsd_amplitudes"))
from torch.utils.data import DataLoader, Subset  # noqa: E402
from pretrain.dataset import CCSDAmplitudeDataset  # noqa: E402
from pretrain.model import ModelConfig, PretrainingModel  # noqa: E402

OUT = Path(__file__).resolve().parent
CACHE = OUT / "art" / "_val_metrics.npz"
BLUE, TEAL, GREY = "#1f6fb2", "#0f8f83", "#98a7b3"


def compute():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dataset = CCSDAmplitudeDataset(
        "rhf_dataset", targets_dir="rhf_targets", inv_targets_dir="rhf_inv_targets")
    gen = torch.Generator().manual_seed(42)
    indices = torch.randperm(len(dataset), generator=gen).tolist()
    n_val = int(len(dataset) * 0.1)
    val_indices = indices[len(dataset) - n_val:]
    loader = DataLoader(Subset(dataset, val_indices), batch_size=12, shuffle=False,
                        collate_fn=dataset.collate_fn, num_workers=4)

    cfg = ModelConfig(embed_dim=192, num_layers=6, num_heads=8, n_reps=2,
                      dropout=0.0, predict_invariant=True)
    model = PretrainingModel(cfg).to(device)
    sd = torch.load("checkpoints_invariant/best.pt", map_location=device,
                    weights_only=False)["model_state_dict"]
    model.load_state_dict(sd, strict=False)
    model.eval()

    lam_t, lam_p, cos, names = [], [], [], []
    maps = {}
    with torch.no_grad():
        for batch in loader:
            b = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                 for k, v in batch.items()}
            out = model(b)
            for i in range(b["t2"].shape[0]):
                no, nv = int(b["noccs"][i]), int(b["nvirts"][i])
                m = out["v_pred"][i][:no, b["max_nocc"]:b["max_nocc"] + nv]
                m = m.cpu().numpy().astype(np.float64)
                v0 = b["v_target"][i][:no, :nv].cpu().numpy().astype(np.float64)
                vh = m / max(np.linalg.norm(m), 1e-12)
                c = float((vh * v0).sum())
                cos.append(abs(c))
                names.append(b["names"][i])
                maps[b["names"][i]] = (v0, np.sign(c) * vh)
                lam_t.append(float(b["lam_target"][i]))
                lam_p.append(float(out["lam_pred"][i]))

    # the median-agreement molecule: a representative case, not the best one
    cos_a = np.array(cos)
    pick = names[int(np.argsort(cos_a)[len(cos_a) // 2])]
    v0, vh = maps[pick]
    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE, lam_t=np.array(lam_t), lam_p=np.array(lam_p), cos=cos_a,
             pick=pick, map_exact=v0, map_pred=vh)
    print(f"cached {len(cos)} molecules -> {CACHE}  (median case: {pick})")


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c9d3db")
    ax.tick_params(colors="#66757f", labelsize=15, length=4)


def render():
    d = np.load(CACHE)
    lam_t, lam_p, cos = np.abs(d["lam_t"]), np.abs(d["lam_p"]), d["cos"]
    r2 = 1 - ((lam_p - lam_t) ** 2).sum() / ((lam_t - lam_t.mean()) ** 2).sum()

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    lo, hi = 0.10, min(np.percentile(lam_t, 99.7), 0.75)
    ax.plot([lo, hi], [lo, hi], color=GREY, lw=2.0, ls=(0, (5, 4)), zorder=1)
    ax.scatter(lam_t, lam_p, s=17, color=TEAL, alpha=0.24,
               edgecolors="none", zorder=2)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("exact calculation", fontsize=19, color="#10344d")
    ax.set_ylabel("AI prediction", fontsize=19, color="#10344d")
    ax.set_xticks([0.2, 0.4, 0.6]); ax.set_yticks([0.2, 0.4, 0.6])
    ax.text(0.045, 0.955, f"R² = {r2:.2f}", transform=ax.transAxes,
            fontsize=27, fontweight="bold", color=TEAL, va="top")
    ax.set_aspect("equal")
    style(ax)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUT / "art" / "parity.png", dpi=200, transparent=True)
    plt.close(fig)

    print(f"n={len(cos)}  R²={r2:.3f}  median cos={np.median(cos):.3f}  "
          f">=0.9: {(cos >= 0.9).mean()*100:.1f}%")


if __name__ == "__main__":
    if not CACHE.exists():
        compute()
    render()
