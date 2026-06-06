#!/usr/bin/env python3
"""Live-updating train/val loss curves from the training log.

Reads the JSON-lines log produced by train.py and plots train vs val
loss curves, refreshing periodically.

Usage:
    python pretrain/plot_curves.py [--log-dir checkpoints] [--interval 30]
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_log(log_path: Path):
    train_steps, train_losses, train_j, train_u = [], [], [], []
    val_epochs, val_losses, val_j, val_u = [], [], [], []

    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue

            if entry.get("type") == "train_step":
                train_steps.append(entry["step"])
                train_losses.append(entry["loss"])
                train_j.append(entry.get("j_loss", 0))
                train_u.append(entry.get("u_loss", 0))
            elif entry.get("type") == "val_epoch":
                val_epochs.append(entry["epoch"])
                val_losses.append(entry["val_loss"])
                val_j.append(entry.get("val_j_loss", 0))
                val_u.append(entry.get("val_u_loss", 0))

    return {
        "train_steps": np.array(train_steps),
        "train_losses": np.array(train_losses),
        "train_j": np.array(train_j),
        "train_u": np.array(train_u),
        "val_epochs": np.array(val_epochs),
        "val_losses": np.array(val_losses),
        "val_j": np.array(val_j),
        "val_u": np.array(val_u),
    }


def smooth(values, window=50):
    if len(values) < window:
        return values
    kernel = np.ones(window) / window
    return np.convolve(values, kernel, mode="valid")


def plot(data, output_path: Path):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    if len(data["train_steps"]) > 0:
        ax = axes[0]
        ax.set_title("Total Loss")
        ax.plot(data["train_steps"], data["train_losses"],
                alpha=0.2, color="C0", label="train (raw)")
        if len(data["train_losses"]) > 50:
            smoothed = smooth(data["train_losses"])
            ax.plot(data["train_steps"][49:], smoothed,
                    color="C0", linewidth=2, label="train (smooth)")
        if len(data["val_epochs"]) > 0:
            ax.plot(data["val_epochs"] * (len(data["train_steps"]) / max(data["val_epochs"].max(), 1)),
                    data["val_losses"], "o-", color="C1", label="val", markersize=4)
        ax.set_xlabel("Step")
        ax.set_ylabel("Loss")
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        ax.set_title("J Loss")
        ax.plot(data["train_steps"], data["train_j"],
                alpha=0.2, color="C0")
        if len(data["train_j"]) > 50:
            ax.plot(data["train_steps"][49:], smooth(data["train_j"]),
                    color="C0", linewidth=2, label="train")
        if len(data["val_j"]) > 0:
            ax.plot(data["val_epochs"] * (len(data["train_steps"]) / max(data["val_epochs"].max(), 1)),
                    data["val_j"], "o-", color="C1", label="val", markersize=4)
        ax.set_xlabel("Step")
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[2]
        ax.set_title("U Loss")
        ax.plot(data["train_steps"], data["train_u"],
                alpha=0.2, color="C0")
        if len(data["train_u"]) > 50:
            ax.plot(data["train_steps"][49:], smooth(data["train_u"]),
                    color="C0", linewidth=2, label="train")
        if len(data["val_u"]) > 0:
            ax.plot(data["val_epochs"] * (len(data["train_steps"]) / max(data["val_epochs"].max(), 1)),
                    data["val_u"], "o-", color="C1", label="val", markersize=4)
        ax.set_xlabel("Step")
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"[{time.strftime('%H:%M:%S')}] Plot saved to {output_path} "
          f"({len(data['train_steps'])} train steps, {len(data['val_epochs'])} val epochs)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--log-dir", default="checkpoints",
                        help="Directory containing train_log.jsonl")
    parser.add_argument("--interval", type=int, default=30,
                        help="Refresh interval in seconds (0 = one-shot)")
    parser.add_argument("--output", default="train_val_curves.png")
    args = parser.parse_args()

    log_path = Path(args.log_dir) / "train_log.jsonl"
    output_path = Path(args.output)

    if not log_path.exists():
        print(f"Waiting for {log_path} to appear...")
        while not log_path.exists():
            time.sleep(5)

    while True:
        data = read_log(log_path)
        if len(data["train_steps"]) > 0:
            plot(data, output_path)
        else:
            print(f"[{time.strftime('%H:%M:%S')}] No data yet...")

        if args.interval <= 0:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
