import sys
import torch
from torch.utils.data import DataLoader, Subset
from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig
from pretrain.train import ReconstructionLoss, BucketBatchSampler

data_dir = sys.argv[1] if len(sys.argv) > 1 else 'rhf_dataset_n29'
ds = CCSDAmplitudeDataset(data_dir, n_reps=2, targets_dir='rhf_targets')

# mirror train.py's split exactly
n_val = int(len(ds) * 0.1)
n_train = len(ds) - n_val
g = torch.Generator().manual_seed(42)
indices = torch.randperm(len(ds), generator=g).tolist()
train_indices, val_indices = indices[:n_train], indices[n_train:]
train_sampler = BucketBatchSampler(ds, 8, shuffle=True, seed=42)
train_norbs = {i: ds.get_norb(i) for i in train_indices}
train_sampler.sorted_indices = sorted(train_indices, key=lambda i: train_norbs[i])
train_loader = DataLoader(ds, batch_sampler=train_sampler, collate_fn=ds.collate_fn, num_workers=0)
val_loader = DataLoader(Subset(ds, val_indices), batch_size=8, shuffle=False,
                        collate_fn=ds.collate_fn, num_workers=0)

dev = 'cuda'
m = PretrainingModel(ModelConfig(embed_dim=192, num_layers=6, num_heads=8, n_reps=2)).to(dev)
crit = ReconstructionLoss(n_reps=2, max_norb=96).to(dev)
opt = torch.optim.AdamW(m.parameters(), lr=3e-3)

print(f"=== iterating VAL loader ({len(val_indices)} mols) ===", flush=True)
for bi, batch in enumerate(val_loader):
    nb = {k: (v.to(dev) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
    try:
        opt.zero_grad()
        out = m(nb)
        r, w = crit(out['kappa_real_pred'], out['kappa_imag_pred'], out['j_pred'],
                    nb['t2'], nb['noccs'], nb['nvirts'], nb['max_nocc'],
                    nb.get('z_target'), nb.get('kappa_real_target'), nb.get('kappa_imag_target'))
        (r + w).backward()
        opt.step()
        torch.cuda.synchronize()
    except Exception as e:
        print('CRASH batch', bi, 'noccs', nb['noccs'].tolist(), 'nvirts', nb['nvirts'].tolist(),
              'max_nocc', nb['max_nocc'], 'max_nvirt', nb['max_nvirt'], flush=True)
        import traceback
        traceback.print_exc()
        break
    if bi % 30 == 0:
        print('val batch', bi, 'ok recon', float(r), 'ws', float(w), flush=True)
print('val done; === iterating TRAIN loader ===', flush=True)

for bi, batch in enumerate(train_loader):
    nb = {k: (v.to(dev) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
    try:
        opt.zero_grad()
        out = m(nb)
        r, w = crit(out['kappa_real_pred'], out['kappa_imag_pred'], out['j_pred'],
                    nb['t2'], nb['noccs'], nb['nvirts'], nb['max_nocc'],
                    nb.get('z_target'), nb.get('kappa_real_target'), nb.get('kappa_imag_target'))
        (r + w).backward()
        opt.step()
        torch.cuda.synchronize()
    except Exception as e:
        print('TRAIN CRASH batch', bi, 'noccs', nb['noccs'].tolist(), 'nvirts', nb['nvirts'].tolist(),
              'norbs', nb['norbs'].tolist(), 'names', nb['names'], flush=True)
        import traceback
        traceback.print_exc()
        break
    if bi % 50 == 0:
        print('train batch', bi, 'ok', float(r), flush=True)
print('done')
