#!/usr/bin/env python3
"""Find which batches produce NaN."""
import torch
from torch.utils.data import DataLoader
from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig

ds = CCSDAmplitudeDataset("jobs")
loader = DataLoader(ds, batch_size=4, collate_fn=ds.collate_fn, num_workers=0)

config = ModelConfig(embed_dim=128, num_layers=1, num_heads=4, n_reps=2)
model = PretrainingModel(config)
model.eval()

for i, batch in enumerate(loader):
    with torch.no_grad():
        t2 = batch["t2"]
        if torch.isnan(t2).any() or torch.isinf(t2).any():
            print(f"Batch {i}: t2 has nan/inf! names={batch['names']}")
            continue

        try:
            out = model(batch)
        except Exception as e:
            print(f"Batch {i}: ERROR {e}, names={batch['names']}")
            continue

        has_nan = torch.isnan(out["j_pred"]).any() or torch.isnan(out["u_real_pred"]).any()
        if has_nan:
            norbs = batch["norbs"]
            nocc = batch["nocc"]
            nvirt = batch["nvirt"]
            print(f"Batch {i}: NaN! norbs={norbs.tolist()}, nocc={nocc}, nvirt={nvirt}, "
                  f"t2_shape={t2.shape}, names={batch['names']}")
            break

    if i >= 100:
        print(f"Checked {i+1} batches, no NaN found in first 100")
        break
    if i % 20 == 0:
        print(f"Batch {i}: OK (norbs={batch['norbs'].tolist()})")
