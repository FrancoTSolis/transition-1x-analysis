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


def build_j_sparsity_mask(norb: int, device: torch.device) -> torch.Tensor:
    """Build the sparsity mask for J matrices.

    J_aa has tridiagonal structure: pairs (p, p+1) for p in [0, norb-1).
    J_ab has diagonal structure: pairs (p, p) for p in [0, norb).

    Returns a (2, norb, norb) boolean mask where dim 0 indexes [aa, ab].
    """
    mask = torch.zeros(2, norb, norb, dtype=torch.bool, device=device)
    # J_aa: tridiagonal (nearest-neighbor)
    for p in range(norb - 1):
        mask[0, p, p + 1] = True
        mask[0, p + 1, p] = True
    # J_ab: diagonal
    for p in range(norb):
        mask[1, p, p] = True
    return mask


class LUCJLoss(nn.Module):
    """Combined loss for J (diag Coulomb) and U (orbital rotation) prediction.

    L = (1/K) * sum_k [ sum_mu ||J_hat - J*||_F^2 (masked)
                       + alpha * sum_mu (||U_re_hat - U_re*||_F^2
                                       + ||U_im_hat - U_im*||_F^2) ]
    """

    def __init__(self, alpha: float = 1.0):
        super().__init__()
        self.alpha = alpha
        self._j_mask_cache: dict[int, torch.Tensor] = {}

    def _get_j_mask(self, norb: int, device: torch.device) -> torch.Tensor:
        if norb not in self._j_mask_cache:
            self._j_mask_cache[norb] = build_j_sparsity_mask(norb, device)
        cached = self._j_mask_cache[norb]
        if cached.device != device:
            cached = cached.to(device)
            self._j_mask_cache[norb] = cached
        return cached

    def forward(
        self,
        j_pred: torch.Tensor,
        j_target: torch.Tensor,
        kappa_real_pred: torch.Tensor,
        kappa_real_target: torch.Tensor,
        kappa_imag_pred: torch.Tensor,
        kappa_imag_target: torch.Tensor,
        norbs: torch.Tensor,
        max_norb: int,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Compute loss on J and kappa with masking for variable-size molecules.

        Returns:
            (total_loss, j_loss, kappa_loss)
        """
        batch_size = j_pred.shape[0]
        device = j_pred.device

        j_loss_total = torch.tensor(0.0, device=device)
        kappa_loss_total = torch.tensor(0.0, device=device)

        for i in range(batch_size):
            n = norbs[i].item()

            j_mask = self._get_j_mask(n, device)
            j_diff = j_pred[i, :, :, :n, :n] - j_target[i, :, :, :n, :n]
            masked_diff = j_diff * j_mask.unsqueeze(0).float()
            j_loss_total = j_loss_total + masked_diff.pow(2).sum()

            kr_diff = kappa_real_pred[i, :, :n, :n] - kappa_real_target[i, :, :n, :n]
            ki_diff = kappa_imag_pred[i, :, :n, :n] - kappa_imag_target[i, :, :n, :n]
            kappa_loss_total = kappa_loss_total + kr_diff.pow(2).sum()
            kappa_loss_total = kappa_loss_total + ki_diff.pow(2).sum()

        j_loss = j_loss_total / batch_size
        kappa_loss = kappa_loss_total / batch_size
        total_loss = j_loss + self.alpha * kappa_loss

        return total_loss, j_loss, kappa_loss


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


def train_one_epoch(
    model: nn.Module,
    loader,
    criterion: LUCJLoss,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LambdaLR,
    scaler: torch.amp.GradScaler | None,
    device: torch.device,
    max_grad_norm: float,
    global_step: int,
    log_interval: int,
    json_log_file=None,
) -> tuple[float, int]:
    """Run one training epoch. Returns (avg_loss, updated_global_step)."""
    model.train()
    total_loss = 0.0
    total_j_loss = 0.0
    total_kappa_loss = 0.0
    num_batches = 0

    for batch in loader:
        batch = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                 for k, v in batch.items()}

        optimizer.zero_grad(set_to_none=True)

        use_amp = scaler is not None
        with torch.amp.autocast("cuda", enabled=use_amp):
            outputs = model(batch)
            loss, j_loss, kappa_loss = criterion(
                j_pred=outputs["j_pred"],
                j_target=batch["j_target"],
                kappa_real_pred=outputs["kappa_real_pred"],
                kappa_real_target=batch["kappa_real_target"],
                kappa_imag_pred=outputs["kappa_imag_pred"],
                kappa_imag_target=batch["kappa_imag_target"],
                norbs=batch["norbs"],
                max_norb=outputs["j_pred"].shape[-1],
            )

        if use_amp:
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), max_grad_norm)
            scaler.step(optimizer)
            scaler.update()
        else:
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_grad_norm)
            optimizer.step()

        scheduler.step()
        global_step += 1

        total_loss += loss.item()
        total_j_loss += j_loss.item()
        total_kappa_loss += kappa_loss.item()
        num_batches += 1

        if global_step % log_interval == 0:
            lr = scheduler.get_last_lr()[0]
            logger.info(
                f"step={global_step:6d}  loss={loss.item():.6f}  "
                f"j_loss={j_loss.item():.6f}  kappa_loss={kappa_loss.item():.6f}  "
                f"lr={lr:.2e}"
            )
            if json_log_file is not None:
                json_log_file.write(json_mod.dumps({
                    "type": "train_step", "step": global_step,
                    "loss": loss.item(), "j_loss": j_loss.item(),
                    "kappa_loss": kappa_loss.item(), "lr": lr,
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
) -> tuple[float, float, float]:
    """Run validation. Returns (avg_total_loss, avg_j_loss, avg_kappa_loss)."""
    model.eval()
    total_loss = 0.0
    total_j_loss = 0.0
    total_kappa_loss = 0.0
    num_batches = 0

    for batch in loader:
        batch = {k: v.to(device) if isinstance(v, torch.Tensor) else v
                 for k, v in batch.items()}

        outputs = model(batch)
        loss, j_loss, kappa_loss = criterion(
            j_pred=outputs["j_pred"],
            j_target=batch["j_target"],
            kappa_real_pred=outputs["kappa_real_pred"],
            kappa_real_target=batch["kappa_real_target"],
            kappa_imag_pred=outputs["kappa_imag_pred"],
            kappa_imag_target=batch["kappa_imag_target"],
            norbs=batch["norbs"],
            max_norb=outputs["j_pred"].shape[-1],
        )

        total_loss += loss.item()
        total_j_loss += j_loss.item()
        total_kappa_loss += kappa_loss.item()
        num_batches += 1

    n = max(num_batches, 1)
    return total_loss / n, total_j_loss / n, total_kappa_loss / n


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
    parser.add_argument("--jobs-dir", type=str, required=True,
                        help="Directory containing job subdirectories with LUCJ data")
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
    parser.add_argument("--alpha", type=float, default=1.0,
                        help="Weight for U (orbital rotation) loss term")
    parser.add_argument("--val-split", type=float, default=0.1,
                        help="Fraction of data for validation")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    parser.add_argument("--checkpoint-dir", type=str, default="checkpoints",
                        help="Directory for saving checkpoints")
    parser.add_argument("--resume", type=str, default=None,
                        help="Path to checkpoint to resume from")
    parser.add_argument("--log-interval", type=int, default=50,
                        help="Log training metrics every N steps")
    parser.add_argument("--num-workers", type=int, default=4,
                        help="DataLoader worker count")
    parser.add_argument("--no-amp", action="store_true",
                        help="Disable automatic mixed precision")

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

    use_amp = torch.cuda.is_available() and not args.no_amp
    if use_amp:
        logger.info("Mixed precision training enabled (torch.amp)")

    # -----------------------------------------------------------------------
    # Data
    # -----------------------------------------------------------------------
    jobs_dir = Path(args.jobs_dir)
    if not jobs_dir.is_dir():
        logger.error(f"Jobs directory not found: {jobs_dir}")
        sys.exit(1)

    dataset = CCSDAmplitudeDataset(jobs_dir)
    logger.info(f"Dataset size: {len(dataset)} molecules")

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

    num_training_steps = args.epochs * len(train_loader)
    scheduler = get_cosine_schedule_with_warmup(
        optimizer, args.warmup_steps, num_training_steps
    )

    scaler = torch.amp.GradScaler("cuda") if use_amp else None
    criterion = LUCJLoss(alpha=args.alpha)

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
                scaler=scaler,
                device=device,
                max_grad_norm=args.max_grad_norm,
                global_step=global_step,
                log_interval=args.log_interval,
                json_log_file=json_log_file,
            )
            epoch_time = time.time() - epoch_start

            # Validation
            val_loss, val_j_loss, val_kappa_loss = validate(
                model, val_loader, criterion, device
            )

            logger.info(
                f"Epoch {epoch+1}/{args.epochs}  "
                f"train_loss={train_loss:.6f}  "
                f"val_loss={val_loss:.6f}  "
                f"val_j={val_j_loss:.6f}  val_kappa={val_kappa_loss:.6f}  "
                f"time={epoch_time:.1f}s"
            )
            json_log_file.write(json_mod.dumps({
                "type": "val_epoch", "epoch": epoch + 1,
                "train_loss": train_loss, "val_loss": val_loss,
                "val_j_loss": val_j_loss, "val_kappa_loss": val_kappa_loss,
                "time": epoch_time,
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
        save_checkpoint(
            checkpoint_dir / "last.pt",
            model, optimizer, scheduler, scaler,
            epoch + 1 if "epoch" in dir() else start_epoch,
            global_step, best_val_loss,
        )
        json_log_file.close()
        logger.info(f"Training finished. Best val loss: {best_val_loss:.6f}")


if __name__ == "__main__":
    main()
