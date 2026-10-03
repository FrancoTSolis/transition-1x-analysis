"""n29 group data for the optimize=True pretraining experiments.

Loads, for the 1410 same-shape molecules (nocc=16, nvirt=13, norb=29):
  t2                      (N, 16, 16, 13, 13)
  canonical exact-DF init  U0 (N, 2, 29, 29) complex, Z0_full (N, 2, 29, 29), znorm_full (N,)
  optimize=True labels     U_opt, Z_opt, resid for a config (square_reg0.005 by default)
from rhf_dataset/ and rhf_targets_compressed/ (stored init = the one the labels were
generated from; it reproduces to 1e-13 when recomputed on another machine).

The fixed split (split_n29_{train,val}.txt) reproduces the split of every earlier
n29 run of pretrain/train.py (seed 42 randperm over the 30205-molecule index, filtered
to the n29 names, last 10% = 141 val molecules), so numbers are comparable to the
0.787 compressed_recon baseline.

Model batches: make_batch(idx) builds the dict PretrainingModel.forward expects
(positional orbital encoding; no padding since all shapes are equal).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[2]          # ccsd_amplitudes/
DATA = ROOT / "rhf_dataset"
LABELS = ROOT / "rhf_targets_compressed"
NAMES = ROOT / "gauge_study" / "names_n29_16_13.txt"
SPLIT_DIR = ROOT / "pretrain" / "opt_true"
NOCC, NVIRT, NORB = 16, 13, 29


def names_n29() -> list[str]:
    return [ln.strip() for ln in open(NAMES) if ln.strip()]


def make_split(write: bool = True):
    """Reproduce train.py's n29 split (compressed_recon runs) and write it to files."""
    index = json.load(open(DATA / "_index.json"))
    all_names = sorted(index.keys())          # CCSDAmplitudeDataset order (no filters)
    gen = torch.Generator().manual_seed(42)
    perm = torch.randperm(len(all_names), generator=gen).tolist()
    allowed = set(names_n29())
    idx = [i for i in perm if all_names[i] in allowed]
    n_val = int(len(idx) * 0.1)
    train = [all_names[i] for i in idx[: len(idx) - n_val]]
    val = [all_names[i] for i in idx[len(idx) - n_val:]]
    if write:
        (SPLIT_DIR / "split_n29_train.txt").write_text("\n".join(train) + "\n")
        (SPLIT_DIR / "split_n29_val.txt").write_text("\n".join(val) + "\n")
    return train, val


def load_split():
    tr = SPLIT_DIR / "split_n29_train.txt"
    va = SPLIT_DIR / "split_n29_val.txt"
    if not (tr.exists() and va.exists()):
        return make_split(True)
    return ([ln.strip() for ln in open(tr) if ln.strip()],
            [ln.strip() for ln in open(va) if ln.strip()])


class N29:
    """All n29 tensors on one device. Index order: self.names."""

    def __init__(self, names: list[str] | None = None, config: str | None = "square_reg0.005",
                 device="cpu", dtype=torch.float32):
        self.names = names or names_n29()
        self.config = config
        cdt = torch.complex64 if dtype == torch.float32 else torch.complex128
        t2, U0, Z0, zf, lam0, v0, gap = [], [], [], [], [], [], []
        Uo, Zo, ro = [], [], []
        for n in self.names:
            t2.append(np.load(DATA / f"{n}.npz")["t2"].astype(np.float64))
            s = np.load(LABELS / "init" / f"{n}.npz")
            U0.append(s["U_re"] + 1j * s["U_im"])
            Z0.append(s["Z"])
            zf.append(float(s["znorm_full"]))
            lam0.append(float(s["lam0"]))
            v0.append(s["v0"])
            gap.append(float(s["gap"]))
            if config:
                c = np.load(LABELS / config / f"{n}.npz")
                Uo.append(c["U_re"] + 1j * c["U_im"])
                Zo.append(c["Z"])
                ro.append(float(c["resid"]))
        T = lambda a, d: torch.as_tensor(np.asarray(a), device=device).to(d)  # noqa: E731
        self.t2 = T(t2, dtype)
        self.U0 = T(U0, cdt)
        self.Z0_full = T(Z0, dtype)
        self.znorm_full = T(zf, dtype)
        self.lam0 = T(lam0, dtype)
        self.v0 = T(v0, dtype)
        self.gap = T(gap, dtype)
        if config:
            self.U_opt = T(Uo, cdt)
            self.Z_opt = T(Zo, dtype)
            self.resid_opt = T(ro, dtype)
        self.device = device
        self.N = len(self.names)
        self.pos = {n: i for i, n in enumerate(self.names)}

    def indices(self, names: list[str]) -> torch.Tensor:
        return torch.tensor([self.pos[n] for n in names], device=self.device)

    def make_batch(self, idx: torch.Tensor, n_reps: int = 2) -> dict:
        """Batch dict for pretrain.model.PretrainingModel (positional mode)."""
        B = len(idx)
        dev = self.t2.device
        return {
            "t2": self.t2[idx].float(),
            "noccs": torch.full((B,), NOCC, device=dev, dtype=torch.long),
            "nvirts": torch.full((B,), NVIRT, device=dev, dtype=torch.long),
            "norbs": torch.full((B,), NORB, device=dev, dtype=torch.long),
            "nocc": [NOCC] * B,
            "nvirt": [NVIRT] * B,
            "n_reps": [n_reps] * B,
            "max_nocc": NOCC,
            "max_nvirt": NVIRT,
        }
