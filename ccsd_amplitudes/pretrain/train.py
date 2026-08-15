#!/usr/bin/env python3
"""Training script for the LUCJ pretraining model.

Trains a model to predict LUCJ diag_coulomb_mats (J) and orbital_rotations (U)
from CCSD amplitudes. Handles variable-size molecules via padding and masking.

Usage:
    python pretrain/train.py --jobs-dir jobs --epochs 200 --batch-size 32
    python pretrain/train.py --jobs-dir jobs --resume checkpoints/best.pt
"""

from __future__ import annotations

import argparse
import json as json_mod
import logging
import math
import signal
import sys
import time
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

from pretrain.dataset import CCSDAmplitudeDataset
from pretrain.model import PretrainingModel, ModelConfig

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------


class LUCJLoss(nn.Module):
    """LUCJ pretraining loss with a selectable U/kappa objective.

    The model predicts, per repetition, an anti-Hermitian generator
    kappa = kappa_real + i*kappa_imag (=> U = expm(kappa)) and a symmetric
    diagonal-Coulomb matrix Z (= j_pred[:, :, 0]). The diagonal-Coulomb J/Z is
    directly learnable from t2 (R^2 ~ 0.79); the orbital rotation U/kappa is
    gauge-scrambled across diverse molecules and is therefore a weak head
    (acceptable for a pretraining stage).

    mode='supervised' (default, robust):
        L = ||Z - Z*||^2 + kappa_weight * (||kr - kr*||^2 + ||ki - ki*||^2)
        Direct regression to the ffsim double-factorization targets. No saddle.
        Inference: U = expm(kappa_pred).

    mode='reconstruction' (gauge-invariant, experimental):
        L = ||t2_rec - t2|| / ||t2||  +  z_anchor_weight * ||Z - Z*||^2
        where t2_rec = i * sum_k,pq Z U U* U U* and U = U_base @ expm(kappa) with
        a FIXED complex base rotation U_base (breaks the (U=I, Z=0) saddle; it
        must be applied at inference too, so it is saved in the checkpoint). The
        Z anchor prevents the Z->0 collapse. This objective is invariant to U's
        gauge, but on heterogeneous data it escapes the saddle only slowly.

    mode='invariant' (deck_v3 retarget — the gauge fix):
        Do not regress kappa at all. Predict the gauge-free sufficient statistic
        of the exact-DF label: lam0 (scalar) and v0 (occ x virt matrix, unit,
        defined up to sign), from which (U, Z) is rebuilt at inference with the
        same deterministic quadrature+eigh construction ffsim uses
        (gauge_study.common.exact_df_t2 / build_from_v).
        L = ||Z - Z*||^2                       (unchanged J/Z head)
          + min(||m - v0||^2, ||m + v0||^2)    (sign is the only residual gauge)
          + lam_weight * (lam - lam0)^2
        Logging channels in this mode: kappa slot = v0 loss, recon slot = lam0
        loss (names kept for compatibility with existing plotting scripts).

    Variable-size molecules: per sample we gather the active orbital slots
    [0..nocc) ∪ [max_nocc..max_nocc+nvirt). ffsim targets are stored in native
    occ-then-virt order, which matches the gathered prediction order.

    Returns (total, recon_val, z_val, kappa_val); the last three are detached
    component values for logging.
    """

    def __init__(self, mode: str = "supervised", kappa_weight: float = 1.0,
                 z_anchor_weight: float = 1.0, recon_relative: bool = True,
                 z_reg: str = "anchor", lam_weight: float = 1.0,
                 n_reps: int = 2, max_norb: int = 96, base_scale: float = 0.3):
        super().__init__()
        assert mode in ("supervised", "reconstruction", "invariant")
        assert z_reg in ("anchor", "floor", "none")
        self.mode = mode
        self.kappa_weight = kappa_weight
        self.lam_weight = lam_weight
        self.z_anchor_weight = z_anchor_weight
        self.z_reg = z_reg
        self.recon_relative = recon_relative
        # Fixed complex base rotation for reconstruction mode (see class docstring).
        g = torch.Generator().manual_seed(0)
        re = base_scale * torch.randn(n_reps, max_norb, max_norb, generator=g)
        im = base_scale * torch.randn(n_reps, max_norb, max_norb, generator=g)
        self.register_buffer("base_re", re - re.transpose(1, 2))   # real antisym
        self.register_buffer("base_im", im + im.transpose(1, 2))   # real sym

    def forward(
        self,
        kappa_real_pred: torch.Tensor,   # (B, reps, N, N) anti-symmetric
        kappa_imag_pred: torch.Tensor,   # (B, reps, N, N) symmetric
        j_pred: torch.Tensor,            # (B, reps, 2, N, N) symmetric (Z = block 0)
        t2: torch.Tensor,                # (B, max_nocc, max_nocc, max_nvirt, max_nvirt)
        noccs: torch.Tensor,
        nvirts: torch.Tensor,
        max_nocc: int,
        z_target: torch.Tensor | None = None,    # (B, reps, N, N) ffsim DF Z
        kr_target: torch.Tensor | None = None,   # (B, reps, N, N) ffsim Re log U
        ki_target: torch.Tensor | None = None,   # (B, reps, N, N) ffsim Im log U
        v_pred: torch.Tensor | None = None,      # (B, N, N) invariant-mode v0 head
        lam_pred: torch.Tensor | None = None,    # (B,) invariant-mode lam0 head
        v_target: torch.Tensor | None = None,    # (B, max_nocc, max_nvirt) unit v0
        lam_target: torch.Tensor | None = None,  # (B,)
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        B = kappa_real_pred.shape[0]
        device = kappa_real_pred.device
        total = torch.zeros((), device=device)
        recon_acc = torch.zeros((), device=device)
        z_acc = torch.zeros((), device=device)
        kappa_acc = torch.zeros((), device=device)

        for i in range(B):
            no = int(noccs[i]); nv = int(nvirts[i]); n = no + nv
            idx = torch.cat([
                torch.arange(no, device=device),
                max_nocc + torch.arange(nv, device=device),
            ])
            kr = kappa_real_pred[i][:, idx][:, :, idx]    # (reps, n, n)
            ki = kappa_imag_pred[i][:, idx][:, :, idx]
            Z = j_pred[i, :, 0][:, idx][:, :, idx]

            if self.mode == "supervised":
                zt = z_target[i][:, :n, :n]
                krt = kr_target[i][:, :n, :n]
                kit = ki_target[i][:, :n, :n]
                z_l = (Z - zt).pow(2).mean()
                k_l = (kr - krt).pow(2).mean() + (ki - kit).pow(2).mean()
                total = total + z_l + self.kappa_weight * k_l
                z_acc = z_acc + z_l.detach()
                kappa_acc = kappa_acc + k_l.detach()
            elif self.mode == "invariant":
                zt = z_target[i][:, :n, :n]
                z_l = (Z - zt).pow(2).mean()
                # active occ-virt block of the pair-token map
                m_hat = v_pred[i][:no, max_nocc:max_nocc + nv]
                y = v_target[i][:no, :nv]
                v_l = torch.minimum(
                    (m_hat - y).pow(2).sum(), (m_hat + y).pow(2).sum())
                lam_l = (lam_pred[i] - lam_target[i]).pow(2)
                total = total + z_l + v_l + self.lam_weight * lam_l
                z_acc = z_acc + z_l.detach()
                kappa_acc = kappa_acc + v_l.detach()   # kappa slot = v0 loss
                recon_acc = recon_acc + lam_l.detach()  # recon slot = lam0 loss
            else:  # reconstruction
                U_head = torch.linalg.matrix_exp(torch.complex(kr, ki))
                b_re = self.base_re[:, idx][:, :, idx]
                b_im = self.base_im[:, idx][:, :, idx]
                U_base = torch.linalg.matrix_exp(torch.complex(b_re, b_im))
                U = U_base @ U_head
                full = 1j * torch.einsum(
                    "kpq,kap,kip,kcq,kjq->ijac",
                    Z.to(U.dtype), U, U.conj(), U, U.conj())
                rec = full[:no, :no, no:, no:].real
                tgt = t2[i, :no, :no, :nv, :nv]
                resid = (rec - tgt).norm()
                if self.recon_relative:
                    resid = resid / tgt.norm().clamp_min(1e-8)
                loss_i = resid
                z_l = torch.zeros((), device=device)
                if z_target is not None and self.z_anchor_weight > 0 and self.z_reg != "none":
                    zt = z_target[i][:, :n, :n]
                    if self.z_reg == "anchor":
                        # pin Z to ffsim's Z (imposes ffsim's gauge on Z)
                        z_l = (Z - zt).pow(2).mean()
                    else:  # "floor": gauge-free magnitude floor
                        # keep ||Z|| from collapsing to 0 without pinning direction,
                        # so U and Z can co-adapt to a consistent reconstructing pair
                        zt_norm = zt.flatten(1).norm(dim=1)
                        z_norm = Z.flatten(1).norm(dim=1)
                        z_l = torch.relu(zt_norm - z_norm).pow(2).mean()
                    loss_i = loss_i + self.z_anchor_weight * z_l
                total = total + loss_i
                recon_acc = recon_acc + resid.detach()
                z_acc = z_acc + z_l.detach()

        return total / B, recon_acc / B, z_acc / B, kappa_acc / B


# ---------------------------------------------------------------------------
# Bucket Sampler
# ---------------------------------------------------------------------------


class BucketBatchSampler:
    """Batch sampler that groups samples by norb to minimize padding waste.

    Sorts indices by norb and yields consecutive chunks of batch_size.
    """

    def __init__(self, dataset: CCSDAmplitudeDataset, batch_size: int, shuffle: bool = True,
                 seed: int = 42):
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.seed = seed
        self.epoch = 0

        norbs = [dataset.get_norb(i) for i in range(len(dataset))]
        self.sorted_indices = sorted(range(len(dataset)), key=lambda i: norbs[i])

    def __iter__(self):
        indices = list(self.sorted_indices)

        if self.shuffle:
            rng = torch.Generator()
            rng.manual_seed(self.seed + self.epoch)
            # Shuffle within buckets of similar size to add randomness
            # while keeping similar norbs together
            bucket_size = self.batch_size * 4
            for start in range(0, len(indices), bucket_size):
                end = min(start + bucket_size, len(indices))
                bucket = indices[start:end]
                perm = torch.randperm(len(bucket), generator=rng).tolist()
                indices[start:end] = [bucket[p] for p in perm]

        batches = []
        for start in range(0, len(indices), self.batch_size):
            batches.append(indices[start:start + self.batch_size])

        if self.shuffle:
            rng = torch.Generator()
            rng.manual_seed(self.seed + self.epoch + 1000)
            perm = torch.randperm(len(batches), generator=rng).tolist()
            batches = [batches[p] for p in perm]

        self.epoch += 1
        return iter(batches)

    def __len__(self):
        return math.ceil(len(self.sorted_indices) / self.batch_size)


# ---------------------------------------------------------------------------
# Scheduler
# ---------------------------------------------------------------------------


def get_cosine_schedule_with_warmup(
    optimizer: torch.optim.Optimizer,
    num_warmup_steps: int,
    num_training_steps: int,
) -> torch.optim.lr_scheduler.LambdaLR:
    """Cosine LR schedule with linear warmup."""

    def lr_lambda(current_step: int) -> float:
        if current_step < num_warmup_steps:
            return current_step / max(1, num_warmup_steps)
        progress = (current_step - num_warmup_steps) / max(
            1, num_training_steps - num_warmup_steps
        )
        return max(0.0, 0.5 * (1.0 + math.cos(math.pi * progress)))

    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------


def _compute_loss(model, batch, criterion, device, return_outputs=False):
    outputs = model(batch)
    total, recon_v, z_v, kappa_v = criterion(
        kappa_real_pred=outputs["kappa_real_pred"],
        kappa_imag_pred=outputs["kappa_imag_pred"],
        j_pred=outputs["j_pred"],
        t2=batch["t2"],
        noccs=batch["noccs"],
        nvirts=batch["nvirts"],
        max_nocc=batch["max_nocc"],
        z_target=batch.get("z_target"),
        kr_target=batch.get("kappa_real_target"),
        ki_target=batch.get("kappa_imag_target"),
        v_pred=outputs.get("v_pred"),
        lam_pred=outputs.get("lam_pred"),
        v_target=batch.get("v_target"),
        lam_target=batch.get("lam_target"),
    )
    if return_outputs:
        return total, recon_v, z_v, kappa_v, outputs
    return total, recon_v, z_v, kappa_v


@torch.no_grad()
def _invariant_batch_metrics(outputs, batch):
    """Per-sample gauge-aware metrics for the invariant mode.

    Returns lists of (cosine, B_rel_err, lam_pred, lam_true) where
    cosine = |<v_hat/||v_hat||, v0>| and
    B_rel_err = ||lam_hat vh vh^T - lam0 v0 v0^T||_F / |lam0| (v's unit).
    """
    cos_l, berr_l, lam_p, lam_t = [], [], [], []
    max_nocc = batch["max_nocc"]
    for i in range(outputs["v_pred"].shape[0]):
        no = int(batch["noccs"][i]); nv = int(batch["nvirts"][i])
        m_hat = outputs["v_pred"][i][:no, max_nocc:max_nocc + nv]
        nrm = m_hat.norm().clamp_min(1e-12)
        vh = m_hat / nrm
        y = batch["v_target"][i][:no, :nv]
        cos = (vh * y).sum().abs().clamp(max=1.0)
        lam_h = outputs["lam_pred"][i]
        lam0 = batch["lam_target"][i]
        # ||B_hat - B||^2 = lam_h^2 + lam0^2 - 2 lam_h lam0 cos^2 (unit v's)
        b2 = (lam_h ** 2 + lam0 ** 2
              - 2.0 * lam_h * lam0 * cos ** 2).clamp_min(0.0)
        berr = b2.sqrt() / lam0.abs().clamp_min(1e-12)
        cos_l.append(cos.item()); berr_l.append(berr.item())
        lam_p.append(lam_h.item()); lam_t.append(lam0.item())
    return cos_l, berr_l, lam_p, lam_t


def train_one_epoch(
    model: nn.Module,
    loader,
    criterion: LUCJLoss,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LambdaLR,
    device: torch.device,
    max_grad_norm: float,
    global_step: int,
    log_interval: int,
    accum_steps: int = 1,
    json_log_file=None,
) -> tuple[float, int]:
    """Run one training epoch with gradient accumulation.

    Returns (avg_total_loss, updated_global_step).
    """
    model.train()
    total_loss = 0.0
    num_batches = 0
    optimizer.zero_grad(set_to_none=True)

    for bi, batch in enumerate(loader):
        batch = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                 for k, v in batch.items()}

        loss, recon_v, z_v, kappa_v = _compute_loss(model, batch, criterion, device)
        (loss / accum_steps).backward()

        total_loss += loss.item()
        num_batches += 1

        if (bi + 1) % accum_steps == 0:
            nn.utils.clip_grad_norm_(model.parameters(), max_grad_norm)
            optimizer.step()
            scheduler.step()
            optimizer.zero_grad(set_to_none=True)
            global_step += 1

            if global_step % log_interval == 0:
                lr = scheduler.get_last_lr()[0]
                logger.info(
                    f"step={global_step:6d}  loss={loss.item():.6f}  "
                    f"z={z_v.item():.6f}  kappa={kappa_v.item():.6f}  "
                    f"recon={recon_v.item():.6f}  lr={lr:.2e}"
                )
                if json_log_file is not None:
                    json_log_file.write(json_mod.dumps({
                        "type": "train_step", "step": global_step,
                        "loss": loss.item(), "z_loss": z_v.item(),
                        "kappa_loss": kappa_v.item(), "recon_loss": recon_v.item(),
                        "lr": lr,
                    }) + "\n")
                    json_log_file.flush()

    avg_loss = total_loss / max(num_batches, 1)
    return avg_loss, global_step


@torch.no_grad()
def validate(
    model: nn.Module,
    loader,
    criterion: LUCJLoss,
    device: torch.device,
) -> tuple[float, float, float, float, dict]:
    """Run validation.

    Returns (total, z, kappa, recon) averages plus a dict of extra
    gauge-aware metrics (populated in invariant mode: median |cos(v0)|,
    median relative B error, lam0 R^2).
    """
    model.eval()
    total_loss = 0.0
    num_batches = 0

    total_z = 0.0
    total_kappa = 0.0
    total_recon = 0.0
    cos_all, berr_all, lam_p_all, lam_t_all = [], [], [], []
    for batch in loader:
        batch = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                 for k, v in batch.items()}
        loss, recon_v, z_v, kappa_v, outputs = _compute_loss(
            model, batch, criterion, device, return_outputs=True)
        total_loss += loss.item()
        total_z += z_v.item()
        total_kappa += kappa_v.item()
        total_recon += recon_v.item()
        num_batches += 1
        if criterion.mode == "invariant" and "v_pred" in outputs:
            c, b, lp, lt = _invariant_batch_metrics(outputs, batch)
            cos_all += c; berr_all += b; lam_p_all += lp; lam_t_all += lt

    n = max(num_batches, 1)
    extras: dict = {}
    if cos_all:
        import numpy as _np
        lam_p = _np.array(lam_p_all); lam_t = _np.array(lam_t_all)
        ss_res = ((lam_p - lam_t) ** 2).sum()
        ss_tot = ((lam_t - lam_t.mean()) ** 2).sum()
        extras = {
            "val_cos_median": float(_np.median(cos_all)),
            "val_cos_mean": float(_np.mean(cos_all)),
            "val_Berr_median": float(_np.median(berr_all)),
            "val_Berr_mean": float(_np.mean(berr_all)),
            "val_lam_r2": float(1.0 - ss_res / max(ss_tot, 1e-12)),
        }
    return total_loss / n, total_z / n, total_kappa / n, total_recon / n, extras


def save_checkpoint(
    path: Path,
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LambdaLR,
    scaler: torch.amp.GradScaler | None,
    epoch: int,
    global_step: int,
    best_val_loss: float,
):
    state = {
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "epoch": epoch,
        "global_step": global_step,
        "best_val_loss": best_val_loss,
    }
    if scaler is not None:
        state["scaler_state_dict"] = scaler.state_dict()
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(state, path)
    logger.info(f"Saved checkpoint to {path}")


def load_checkpoint(
    path: Path,
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LambdaLR,
    scaler: torch.amp.GradScaler | None,
    device: torch.device,
) -> tuple[int, int, float]:
    """Load checkpoint. Returns (epoch, global_step, best_val_loss)."""
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model_state_dict"])
    optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
    scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
    if scaler is not None and "scaler_state_dict" in checkpoint:
        scaler.load_state_dict(checkpoint["scaler_state_dict"])
    epoch = checkpoint["epoch"]
    global_step = checkpoint["global_step"]
    best_val_loss = checkpoint["best_val_loss"]
    logger.info(
        f"Resumed from {path} (epoch={epoch}, step={global_step}, "
        f"best_val={best_val_loss:.6f})"
    )
    return epoch, global_step, best_val_loss


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train the LUCJ pretraining model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=str, default="rhf_dataset",
                        help="Directory of per-molecule RHF-CCSD .npz files "
                             "(from generate_rhf_dataset.py)")
    parser.add_argument("--epochs", type=int, default=200,
                        help="Number of training epochs")
    parser.add_argument("--batch-size", type=int, default=32,
                        help="Training batch size")
    parser.add_argument("--lr", type=float, default=3e-4,
                        help="Peak learning rate")
    parser.add_argument("--weight-decay", type=float, default=1e-5,
                        help="AdamW weight decay")
    parser.add_argument("--warmup-steps", type=int, default=1000,
                        help="Linear warmup steps for LR schedule")
    parser.add_argument("--max-grad-norm", type=float, default=1.0,
                        help="Max gradient norm for clipping")
    parser.add_argument("--accum-steps", type=int, default=1,
                        help="Gradient accumulation steps (effective batch = "
                             "batch-size * accum-steps; cuts gradient noise)")
    parser.add_argument("--targets-dir", type=str, default="rhf_targets",
                        help="Directory of ffsim DF (Z, kappa) targets")
    parser.add_argument("--loss-mode", type=str, default="supervised",
                        choices=["supervised", "reconstruction", "invariant"],
                        help="U/kappa objective: 'supervised' = direct regression "
                             "to ffsim (Z, kappa) targets (robust; J/Z learns well, "
                             "kappa head is weak); 'reconstruction' = gauge-invariant "
                             "t2 reconstruction + Z anchor (experimental); "
                             "'invariant' = deck_v3 retarget: predict (lam0, v0) "
                             "with a sign-invariant loss, rebuild (U, Z) at "
                             "inference via quadratures+eigh (no gauge in the "
                             "objective).")
    parser.add_argument("--kappa-weight", type=float, default=1.0,
                        help="[supervised] weight of the kappa regression term")
    parser.add_argument("--inv-targets-dir", type=str, default="rhf_inv_targets",
                        help="[invariant] directory of (lam0, v0) targets from "
                             "build_invariant_targets.py")
    parser.add_argument("--lam-weight", type=float, default=1.0,
                        help="[invariant] weight of the lam0 regression term")
    parser.add_argument("--z-anchor-weight", type=float, default=1.0,
                        help="[reconstruction] weight of the Z regularization term")
    parser.add_argument("--z-reg", type=str, default="anchor",
                        choices=["anchor", "floor", "none"],
                        help="[reconstruction] Z regularizer: 'anchor' pins Z to "
                             "ffsim's Z (imposes its gauge); 'floor' is a gauge-free "
                             "magnitude floor (prevents Z->0 collapse, lets U,Z "
                             "co-adapt); 'none' = pure reconstruction")
    parser.add_argument("--val-split", type=float, default=0.1,
                        help="Fraction of data for validation")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    parser.add_argument("--checkpoint-dir", type=str, default="checkpoints",
                        help="Directory for saving checkpoints")
    parser.add_argument("--resume", type=str, default=None,
                        help="Path to checkpoint to resume from (restores optimizer/epoch)")
    parser.add_argument("--init-from", type=str, default=None,
                        help="Warm-start: load ONLY model weights from this checkpoint "
                             "(fresh optimizer/epoch). Use to init a reconstruction run "
                             "from a supervised checkpoint so Z starts nonzero/structured "
                             "and the model is out of the Z=0 trap.")
    parser.add_argument("--log-interval", type=int, default=50,
                        help="Log training metrics every N steps")
    parser.add_argument("--num-workers", type=int, default=4,
                        help="DataLoader worker count")

    # Model config overrides
    parser.add_argument("--embed-dim", type=int, default=192,
                        help="Model embedding dimension")
    parser.add_argument("--num-layers", type=int, default=6,
                        help="Number of transformer layers")
    parser.add_argument("--num-heads", type=int, default=8,
                        help="Number of attention heads")
    parser.add_argument("--n-reps", type=int, default=2,
                        help="Number of LUCJ repetitions (layers)")
    parser.add_argument("--dropout", type=float, default=0.1,
                        help="Dropout rate")
    parser.add_argument("--use-hf-energies", action="store_true",
                        help="Use HF orbital energies instead of learned positional embeddings")
    parser.add_argument("--use-mo-coeffs", action="store_true",
                        help="Use atom-aware set pooling of MO coefficients for orbital encoding")
    parser.add_argument("--species-filter", type=str, default=None,
                        choices=["TS", "RP"],
                        help="Filter data by species: TS=transition states only, RP=reactants+products only")
    return parser.parse_args()


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger.info(f"Arguments: {vars(args)}")

    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    # AMP is disabled: the reconstruction loss uses complex matrix_exp, which
    # requires float32 precision.
    scaler = None

    # -----------------------------------------------------------------------
    # Data
    # -----------------------------------------------------------------------
    data_dir = Path(args.data_dir)
    if not data_dir.is_dir():
        logger.error(f"Data directory not found: {data_dir}")
        sys.exit(1)

    dataset = CCSDAmplitudeDataset(
        data_dir, species_filter=args.species_filter, n_reps=args.n_reps,
        targets_dir=args.targets_dir,
        inv_targets_dir=(args.inv_targets_dir
                         if args.loss_mode == "invariant" else None),
    )
    logger.info(f"Dataset size: {len(dataset)} molecules"
                + (f" (with Z targets from {args.targets_dir})" if args.targets_dir else ""))

    n_val = int(len(dataset) * args.val_split)
    n_train = len(dataset) - n_val

    generator = torch.Generator().manual_seed(args.seed)
    indices = torch.randperm(len(dataset), generator=generator).tolist()
    train_indices = indices[:n_train]
    val_indices = indices[n_train:]

    train_dataset = Subset(dataset, train_indices)
    val_dataset = Subset(dataset, val_indices)
    logger.info(f"Train: {n_train}, Val: {n_val}")

    # Bucket sampler for training (groups by norb)
    train_sampler = BucketBatchSampler(
        dataset, args.batch_size, shuffle=True, seed=args.seed
    )
    # Filter sampler indices to only include train set
    train_norbs = {idx: dataset.get_norb(idx) for idx in train_indices}
    train_sampler.sorted_indices = sorted(train_indices, key=lambda i: train_norbs[i])

    train_loader = DataLoader(
        dataset,
        batch_sampler=train_sampler,
        collate_fn=dataset.collate_fn,
        num_workers=args.num_workers,
        pin_memory=torch.cuda.is_available(),
    )

    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        collate_fn=dataset.collate_fn,
        num_workers=args.num_workers,
        pin_memory=torch.cuda.is_available(),
    )

    # -----------------------------------------------------------------------
    # Model
    # -----------------------------------------------------------------------
    model_config = ModelConfig(
        embed_dim=args.embed_dim,
        num_layers=args.num_layers,
        num_heads=args.num_heads,
        n_reps=args.n_reps,
        dropout=args.dropout,
        use_hf_energies=args.use_hf_energies,
        use_mo_coeffs=args.use_mo_coeffs,
        predict_invariant=(args.loss_mode == "invariant"),
    )
    model = PretrainingModel(model_config).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info(f"Model parameters: {n_params:,}")

    # -----------------------------------------------------------------------
    # Optimizer, scheduler, loss
    # -----------------------------------------------------------------------
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay
    )

    num_training_steps = args.epochs * (len(train_loader) // args.accum_steps)
    scheduler = get_cosine_schedule_with_warmup(
        optimizer, args.warmup_steps, num_training_steps
    )

    criterion = LUCJLoss(
        mode=args.loss_mode,
        kappa_weight=args.kappa_weight,
        z_anchor_weight=args.z_anchor_weight,
        z_reg=args.z_reg,
        lam_weight=args.lam_weight,
        n_reps=args.n_reps,
        max_norb=model_config.max_norb,
    ).to(device)
    logger.info(f"Loss mode: {args.loss_mode}"
                + (f" (z_reg={args.z_reg})" if args.loss_mode == "reconstruction" else ""))

    # Warm-start: load only model weights (fresh optimizer/epoch).
    if args.init_from:
        ipath = Path(args.init_from)
        if not ipath.exists():
            logger.error(f"--init-from checkpoint not found: {ipath}")
            sys.exit(1)
        sd = torch.load(ipath, map_location=device, weights_only=False)["model_state_dict"]
        missing, unexpected = model.load_state_dict(sd, strict=False)
        logger.info(f"Warm-started model weights from {ipath} "
                    f"(missing={len(missing)}, unexpected={len(unexpected)})")

    # -----------------------------------------------------------------------
    # Resume
    # -----------------------------------------------------------------------
    start_epoch = 0
    global_step = 0
    best_val_loss = float("inf")

    if args.resume:
        resume_path = Path(args.resume)
        if not resume_path.exists():
            logger.error(f"Checkpoint not found: {resume_path}")
            sys.exit(1)
        start_epoch, global_step, best_val_loss = load_checkpoint(
            resume_path, model, optimizer, scheduler, scaler, device
        )

    # -----------------------------------------------------------------------
    # Training loop
    # -----------------------------------------------------------------------
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    interrupted = False

    def handle_interrupt(signum, frame):
        nonlocal interrupted
        if interrupted:
            logger.warning("Second interrupt received, exiting immediately")
            sys.exit(1)
        interrupted = True
        logger.warning("Interrupt received, finishing current step and saving...")

    signal.signal(signal.SIGINT, handle_interrupt)
    signal.signal(signal.SIGTERM, handle_interrupt)

    logger.info(
        f"Starting training: epochs={args.epochs}, steps/epoch~={len(train_loader)}, "
        f"total_steps~={num_training_steps}"
    )

    json_log_path = checkpoint_dir / "train_log.jsonl"
    json_log_file = open(json_log_path, "a")

    try:
        for epoch in range(start_epoch, args.epochs):
            if interrupted:
                break

            epoch_start = time.time()
            train_loss, global_step = train_one_epoch(
                model=model,
                loader=train_loader,
                criterion=criterion,
                optimizer=optimizer,
                scheduler=scheduler,
                device=device,
                max_grad_norm=args.max_grad_norm,
                global_step=global_step,
                log_interval=args.log_interval,
                accum_steps=args.accum_steps,
                json_log_file=json_log_file,
            )
            epoch_time = time.time() - epoch_start

            # Validation
            val_loss, val_z, val_kappa, val_recon, val_extras = validate(
                model, val_loader, criterion, device)

            extra_str = "".join(
                f"  {k}={v:.4f}" for k, v in sorted(val_extras.items()))
            logger.info(
                f"Epoch {epoch+1}/{args.epochs}  "
                f"train_loss={train_loss:.6f}  val_loss={val_loss:.6f}  "
                f"val_z={val_z:.6f}  val_kappa={val_kappa:.6f}  "
                f"val_recon={val_recon:.6f}{extra_str}  time={epoch_time:.1f}s"
            )
            json_log_file.write(json_mod.dumps({
                "type": "val_epoch", "epoch": epoch + 1,
                "train_loss": train_loss, "val_loss": val_loss,
                "val_z_loss": val_z, "val_kappa_loss": val_kappa,
                "val_recon_loss": val_recon, "time": epoch_time,
                **val_extras,
            }) + "\n")
            json_log_file.flush()

            # Save best model
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                save_checkpoint(
                    checkpoint_dir / "best.pt",
                    model, optimizer, scheduler, scaler,
                    epoch + 1, global_step, best_val_loss,
                )

            # Periodic checkpoint
            if (epoch + 1) % 10 == 0:
                save_checkpoint(
                    checkpoint_dir / f"epoch_{epoch+1:04d}.pt",
                    model, optimizer, scheduler, scaler,
                    epoch + 1, global_step, best_val_loss,
                )

            if interrupted:
                break

    except KeyboardInterrupt:
        logger.warning("KeyboardInterrupt caught")
    finally:
        try:
            save_checkpoint(
                checkpoint_dir / "last.pt",
                model, optimizer, scheduler, scaler,
                epoch + 1 if "epoch" in dir() else start_epoch,
                global_step, best_val_loss,
            )
        except Exception as e:
            logger.error(f"Could not save final checkpoint: {e}")
        json_log_file.close()
        logger.info(f"Training finished. Best val loss: {best_val_loss:.6f}")


if __name__ == "__main__":
    main()
