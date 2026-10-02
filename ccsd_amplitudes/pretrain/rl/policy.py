"""Residual LUCJ policy for energy-based RL fine-tuning.

The input t2 determines the canonical exact-DF init (U_init, Z_init) exactly
(gauge_study.compressed_canonical.canonical_exact_init) -- no prediction is
needed for it.  The value the network adds is a *residual* on top of it:

    U_k  = U_init,k @ expm(dkappa_k)          dkappa anti-Hermitian
    Z_k  = Z_init,k + dZ_k                    dZ real symmetric, LUCJ-masked

The Edge Transformer backbone (pretrain.model.PretrainingModel) produces pair
tokens; ResidualHeads reads dkappa/dZ out of them exactly like the original
kappa/J heads (anti-symmetric / symmetric by construction).  With zero-init
heads the initial policy *is* the canonical exact-DF init.

Gaussian policy for GRPO:  a = mu_theta(t2) + sigma * eps  on the flattened
free parameters (upper triangles), so that
    log pi(a | t2) = -||a - mu||^2 / (2 sigma^2) + const,
    log ratio      = (||a - mu_old||^2 - ||a - mu||^2) / (2 sigma^2),
    KL(pi || pi_ref) = ||mu - mu_ref||^2 / (2 sigma^2).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn
from torch import Tensor

from pretrain.model import ModelConfig, PretrainingModel, ResidualHeads


@dataclass
class PolicyConfig:
    n_reps: int = 2
    kappa_scale: float = 1.0       # tanh bound on dkappa entries
    z_scale: float = 1.0           # tanh bound on dZ entries
    zero_init: bool = True         # start exactly at the canonical init
    connectivity: str = "square"


class ResidualPolicy(nn.Module):
    """Backbone from a (pretrained) PretrainingModel + residual heads."""

    def __init__(self, model_cfg: ModelConfig, pol_cfg: PolicyConfig):
        super().__init__()
        self.base = PretrainingModel(model_cfg)
        self.heads = ResidualHeads(model_cfg.embed_dim, pol_cfg.n_reps,
                                   kappa_scale=pol_cfg.kappa_scale,
                                   z_scale=pol_cfg.z_scale, zero_init=pol_cfg.zero_init)
        self.pol_cfg = pol_cfg

    def load_backbone(self, state_dict: dict, strict: bool = False):
        missing, unexpected = self.base.load_state_dict(state_dict, strict=strict)
        return missing, unexpected

    def backbone_tokens(self, batch) -> tuple[Tensor, int]:
        m = self.base
        t2 = batch["t2"]
        nocc, nvirt = t2.shape[1], t2.shape[3]
        norb = nocc + nvirt
        n_reps = batch["n_reps"][0] if isinstance(batch["n_reps"], list) else batch["n_reps"]
        x = m.tokenizer(t2, nocc, nvirt, n_reps, None,
                        orb_energies=batch.get("orb_energies"),
                        mo_coeffs=batch.get("mo_coeffs"), ao_Z=batch.get("ao_Z"),
                        ao_l=batch.get("ao_l"))
        t2_bias = m._compute_t2_bias(t2, nocc, nvirt, norb)
        norbs = batch.get("norbs")
        mask = m._compute_batch_mask(norbs, norb, x.device) if norbs is not None else None
        for layer in m.layers:
            x = layer(x, mask=mask, t2_bias=t2_bias)
        return m.final_norm(x), norb

    def forward(self, batch) -> dict[str, Tensor]:
        x, norb = self.backbone_tokens(batch)
        return self.heads(x, norb)


# ------------------------------------------------------ flat action space

def free_index(norb: int, z_mask: np.ndarray, n_reps: int = 2):
    """Indices of the free real parameters of (dkappa_re antisym, dkappa_im sym,
    dZ sym-masked) inside the padded (R, N, N) arrays, for one molecule."""
    iu1 = np.triu_indices(norb, 1)
    iu0 = np.triu_indices(norb, 0)
    zr, zc = np.nonzero(np.triu(z_mask))
    return dict(kr=iu1, ki=iu0, dz=(zr, zc))


def gather_flat(out: dict[str, Tensor], i: int, idx_active: Tensor, fidx: dict) -> Tensor:
    """Flatten sample i's free parameters (active slots idx_active) to a vector."""
    kr = out["dkappa_re"][i][:, idx_active][:, :, idx_active]
    ki = out["dkappa_im"][i][:, idx_active][:, :, idx_active]
    dz = out["dz"][i][:, idx_active][:, :, idx_active]
    return torch.cat([kr[:, fidx["kr"][0], fidx["kr"][1]].reshape(-1),
                      ki[:, fidx["ki"][0], fidx["ki"][1]].reshape(-1),
                      dz[:, fidx["dz"][0], fidx["dz"][1]].reshape(-1)])


def unflatten_np(a: np.ndarray, norb: int, fidx: dict, n_reps: int = 2):
    """Inverse of gather_flat for one molecule (numpy)."""
    n1 = n_reps * len(fidx["kr"][0]); n0 = n_reps * len(fidx["ki"][0])
    kr_v = a[:n1].reshape(n_reps, -1); ki_v = a[n1:n1 + n0].reshape(n_reps, -1)
    dz_v = a[n1 + n0:].reshape(n_reps, -1)
    dk = np.zeros((n_reps, norb, norb), dtype=complex)
    dZ = np.zeros((n_reps, norb, norb))
    for r in range(n_reps):
        A = np.zeros((norb, norb)); A[fidx["kr"]] = kr_v[r]; A = A - A.T
        S = np.zeros((norb, norb)); S[fidx["ki"]] = ki_v[r]; S = S + S.T - np.diag(np.diag(S))
        dk[r] = A + 1j * S
        Zr = np.zeros((norb, norb)); Zr[fidx["dz"]] = dz_v[r]
        Zr = Zr + Zr.T - np.diag(np.diag(Zr))
        dZ[r] = Zr
    return dk, dZ


def apply_residual(U_init: np.ndarray, Z_init: np.ndarray, dk: np.ndarray, dZ: np.ndarray):
    """(U_init, Z_init) + residual -> (U, Z)."""
    import scipy.linalg
    U = np.stack([U_init[r] @ scipy.linalg.expm(dk[r]) for r in range(U_init.shape[0])])
    return U, Z_init + dZ
