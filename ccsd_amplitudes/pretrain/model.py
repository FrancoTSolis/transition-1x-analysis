from __future__ import annotations

import math
from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor


@dataclass
class ModelConfig:
    embed_dim: int = 192
    num_layers: int = 6
    num_heads: int = 8
    ffn_multiplier: int = 4
    dropout: float = 0.1
    attention_dropout: float = 0.1
    max_norb: int = 96
    n_reps: int = 2
    use_hf_energies: bool = False
    use_mo_coeffs: bool = False


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------


class DeepSetsEncoder(nn.Module):
    """Encode variable-length flattened slices into fixed-dim vectors."""

    def __init__(self, hidden_dim: int, output_dim: int):
        super().__init__()
        self.phi = nn.Sequential(
            nn.Linear(1, hidden_dim),
            nn.GELU(),
        )
        self.rho = nn.Sequential(
            nn.Linear(hidden_dim, output_dim),
            nn.GELU(),
        )

    def forward(self, x: Tensor) -> Tensor:
        """x: (..., L) variable-length flattened slice values → (..., output_dim)."""
        h = self.phi(x.unsqueeze(-1))  # (..., L, hidden)
        h = h.mean(dim=-2)  # (..., hidden)
        return self.rho(h)  # (..., output_dim)


# ---------------------------------------------------------------------------
# Tokenizer
# ---------------------------------------------------------------------------


class AtomAwareOrbitalEncoder(nn.Module):
    """Encode each MO via atom-aware set pooling over its AO coefficients.

    For each MO p, its coefficient vector C[p, :] has one entry per AO.
    Each AO mu has metadata: nuclear charge Z_mu and angular momentum l_mu.
    We encode each AO contribution as phi(c_{p,mu}, Z_mu, l_mu), then
    mean-pool across AOs to get a fixed-size orbital embedding.
    """

    def __init__(self, output_dim: int, max_Z: int = 54, max_l: int = 4):
        super().__init__()
        self.Z_embed = nn.Embedding(max_Z + 1, 16)
        self.l_embed = nn.Embedding(max_l + 1, 8)
        self.phi = nn.Sequential(
            nn.Linear(1 + 16 + 8, output_dim),
            nn.GELU(),
            nn.Linear(output_dim, output_dim),
        )

    def forward(
        self, mo_coeffs: Tensor, ao_Z: Tensor, ao_l: Tensor, norb: int
    ) -> Tensor:
        B, max_norb, max_nbasis = mo_coeffs.shape
        coeffs = mo_coeffs[:, :norb]  # (B, norb, nbasis)

        Z_feats = self.Z_embed(ao_Z.long())  # (B, nbasis, 16)
        l_feats = self.l_embed(ao_l.long())   # (B, nbasis, 8)
        ao_feats = torch.cat([Z_feats, l_feats], dim=-1)  # (B, nbasis, 24)

        ao_feats_exp = ao_feats.unsqueeze(1).expand(-1, norb, -1, -1)
        coeffs_exp = coeffs.unsqueeze(-1)
        combined = torch.cat([coeffs_exp, ao_feats_exp], dim=-1)  # (B, norb, nbasis, 25)

        encoded = self.phi(combined)  # (B, norb, nbasis, output_dim)
        ao_mask = (ao_Z > 0).unsqueeze(1).unsqueeze(-1)  # (B, 1, nbasis, 1)
        encoded = encoded * ao_mask
        count = ao_mask.sum(dim=2).clamp(min=1)
        pooled = encoded.sum(dim=2) / count  # (B, norb, output_dim)
        return pooled


class Tokenizer(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        d = config.embed_dim
        self.use_hf_energies = config.use_hf_energies
        self.use_mo_coeffs = config.use_mo_coeffs

        # (a) Per-orbital: type embedding (occ=0, virt=1)
        self.type_embed = nn.Embedding(2, d // 4)

        if self.use_mo_coeffs:
            self.mo_encoder = AtomAwareOrbitalEncoder(output_dim=d // 4)
        elif self.use_hf_energies:
            self.energy_proj = nn.Sequential(
                nn.Linear(1, d // 4),
                nn.GELU(),
                nn.Linear(d // 4, d // 4),
            )
        else:
            self.pos_embed = nn.Embedding(config.max_norb, d // 4)

        self.orbital_proj = nn.Sequential(
            nn.Linear(d, d),
            nn.GELU(),
            nn.Linear(d, d),
        )

        # (b) Pair-type embedding: {occ-occ, occ-virt, virt-occ, virt-virt}
        self.pair_type_embed = nn.Embedding(4, d)

        # (c) t2-slice encoder (DeepSets)
        self.slice_encoder = DeepSetsEncoder(hidden_dim=d // 2, output_dim=d)

        # (d) Global features: [nocc, nvirt, n_reps] → d
        self.global_proj = nn.Linear(3, d)

    def forward(
        self,
        t2: Tensor,
        nocc: int,
        nvirt: int,
        n_reps: int,
        pad_norb: int | None = None,
        orb_energies: Tensor | None = None,
        mo_coeffs: Tensor | None = None,
        ao_Z: Tensor | None = None,
        ao_l: Tensor | None = None,
    ) -> Tensor:
        B = t2.shape[0]
        norb = nocc + nvirt
        N = pad_norb if pad_norb is not None else norb
        device = t2.device

        # --- (a) Orbital feature embedding ---
        orb_types = torch.cat([
            torch.zeros(nocc, dtype=torch.long, device=device),
            torch.ones(nvirt, dtype=torch.long, device=device),
        ])

        type_feats = self.type_embed(orb_types)  # (norb, d//4)
        type_feats_b = type_feats.unsqueeze(0).expand(B, -1, -1)

        if self.use_mo_coeffs and mo_coeffs is not None:
            mo_feats = self.mo_encoder(mo_coeffs, ao_Z, ao_l, norb)
            orb_feats = torch.cat([type_feats_b, mo_feats], dim=-1)
        elif self.use_hf_energies and orb_energies is not None:
            ene = orb_energies[:, :norb].unsqueeze(-1)
            ene_feats = self.energy_proj(ene)
            orb_feats = torch.cat([type_feats_b, ene_feats], dim=-1)
        else:
            orb_pos = torch.arange(norb, device=device)
            pos_feats = self.pos_embed(orb_pos)
            orb_feats = torch.cat([type_feats, pos_feats], dim=-1)
            orb_feats = orb_feats.unsqueeze(0).expand(B, -1, -1)

        # Pair (p, q): concat orbital features — orb_feats is (B, norb, d//2)
        p_feats = orb_feats.unsqueeze(2).expand(-1, -1, norb, -1)  # (B, norb, norb, d//2)
        q_feats = orb_feats.unsqueeze(1).expand(-1, norb, -1, -1)  # (B, norb, norb, d//2)
        pair_orb = self.orbital_proj(
            torch.cat([p_feats, q_feats], dim=-1)
        )  # (B, norb, norb, d)

        # --- (b) Pair-type embedding ---
        p_is_virt = orb_types.unsqueeze(1)  # (norb, 1)
        q_is_virt = orb_types.unsqueeze(0)  # (1, norb)
        pair_type_idx = p_is_virt * 2 + q_is_virt  # (norb, norb)
        pair_type_feats = self.pair_type_embed(pair_type_idx)  # (norb, norb, d)
        pair_type_feats = pair_type_feats.unsqueeze(0).expand(B, -1, -1, -1)

        # --- (c) t2-slice features ---
        slice_feats = self._encode_t2_slices(t2, nocc, nvirt, norb)

        # --- (d) Global features ---
        global_vec = torch.tensor(
            [[float(nocc), float(nvirt), float(n_reps)]],
            device=device, dtype=torch.float32,
        ).expand(B, -1)
        global_feats = self.global_proj(global_vec)[:, None, None, :]  # (B, 1, 1, d)

        # --- Combine ---
        tokens = pair_orb + pair_type_feats + slice_feats + global_feats

        # Pad to N if needed
        if N > norb:
            tokens = F.pad(tokens, (0, 0, 0, N - norb, 0, N - norb))

        return tokens

    def _encode_t2_slices(
        self, t2: Tensor, nocc: int, nvirt: int, norb: int
    ) -> Tensor:
        B = t2.shape[0]
        device = t2.device
        d = self.slice_encoder.rho[0].out_features

        result = torch.zeros(B, norb, norb, d, device=device)

        # occ-occ: slice = t2[p, q, :, :] ∈ R^{nvirt×nvirt}
        oo = t2.reshape(B, nocc, nocc, nvirt * nvirt)
        result[:, :nocc, :nocc] = self.slice_encoder(oo)

        # occ-virt: slice = t2[p, :, q-nocc, :] ∈ R^{nocc×nvirt}
        ov = t2.permute(0, 1, 3, 2, 4).reshape(B, nocc, nvirt, nocc * nvirt)
        result[:, :nocc, nocc:norb] = self.slice_encoder(ov)

        # virt-virt: slice = t2[:, :, p-nocc, q-nocc] ∈ R^{nocc×nocc}
        vv = t2.permute(0, 3, 4, 1, 2).reshape(B, nvirt, nvirt, nocc * nocc)
        result[:, nocc:norb, nocc:norb] = self.slice_encoder(vv)

        # virt-occ: slice = t2[:, q, p-nocc, :] ∈ R^{nocc×nvirt}
        vo = t2.permute(0, 3, 2, 1, 4).reshape(B, nvirt, nocc, nocc * nvirt)
        result[:, nocc:norb, :nocc] = self.slice_encoder(vo)

        return result


# ---------------------------------------------------------------------------
# Edge Attention with t2 bias
# ---------------------------------------------------------------------------


class EdgeAttentionWithBias(nn.Module):
    """Edge Transformer attention with additive t2 structural bias.

    Operates on pair tokens of shape (B, N, N, d). Uses the triangular
    attention pattern: score[x, a, y] measures pair (x,a) attending to (a,y).
    """

    def __init__(self, embed_dim: int, num_heads: int, dropout: float):
        super().__init__()
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.d_k = embed_dim // num_heads

        self.left_k_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.right_k_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.left_v_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.right_v_proj = nn.Linear(embed_dim, embed_dim, bias=False)
        self.out_proj = nn.Linear(embed_dim, embed_dim, bias=False)

        self.attn_dropout = nn.Dropout(dropout)
        self.bias_scale = nn.Parameter(torch.ones(1))

    def forward(
        self, x: Tensor, mask: Tensor | None = None, t2_bias: Tensor | None = None
    ) -> Tensor:
        B, N1, N2, _ = x.shape
        H, d_k = self.num_heads, self.d_k

        left_k = self.left_k_proj(x).view(B, N1, N2, H, d_k)
        right_k = self.right_k_proj(x).view(B, N1, N2, H, d_k)
        left_v = self.left_v_proj(x).view(B, N1, N2, H, d_k)
        right_v = self.right_v_proj(x).view(B, N1, N2, H, d_k)

        scores = torch.einsum("bxahd,bayhd->bxayh", left_k, right_k) / math.sqrt(d_k)

        if t2_bias is not None:
            scores = scores + self.bias_scale * t2_bias

        if mask is not None:
            scores = scores.masked_fill(mask.unsqueeze(-1), float("-inf"))

        attn = F.softmax(scores, dim=2)
        attn = attn.nan_to_num(0.0)
        attn = self.attn_dropout(attn)

        val = torch.einsum("bxahd,bayhd->bxayhd", left_v, right_v)
        out = torch.einsum("bxayh,bxayhd->bxyhd", attn, val)
        out = out.reshape(B, N1, N2, self.embed_dim)

        return self.out_proj(out)


# ---------------------------------------------------------------------------
# Transformer layer (pre-norm)
# ---------------------------------------------------------------------------


class TransformerLayer(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        d = config.embed_dim
        self.norm1 = nn.LayerNorm(d)
        self.attn = EdgeAttentionWithBias(d, config.num_heads, config.attention_dropout)
        self.norm2 = nn.LayerNorm(d)
        self.ffn = nn.Sequential(
            nn.Linear(d, d * config.ffn_multiplier),
            nn.GELU(),
            nn.Dropout(config.dropout),
            nn.Linear(d * config.ffn_multiplier, d),
            nn.Dropout(config.dropout),
        )

    def forward(
        self, x: Tensor, mask: Tensor | None = None, t2_bias: Tensor | None = None
    ) -> Tensor:
        h = self.norm1(x)
        x = x + self.attn(h, mask=mask, t2_bias=t2_bias)
        h = self.norm2(x)
        x = x + self.ffn(h)
        return x


# ---------------------------------------------------------------------------
# Decode heads
# ---------------------------------------------------------------------------


class DecodeHeads(nn.Module):
    """Predict LUCJ J and U matrices from final pair embeddings.

    J head: symmetric via g(z_pq) + g(z_qp), then sparsity-masked.
    U head: predict real and imaginary parts of unitary matrices.
    """

    def __init__(self, config: ModelConfig):
        super().__init__()
        d = config.embed_dim
        n_reps = config.n_reps

        # J: one MLP per (rep, spin_type). 2 spin types → 2*n_reps heads.
        n_j_heads = 2 * n_reps
        self.j_mlps = nn.ModuleList([
            nn.Sequential(nn.Linear(d, d // 2), nn.GELU(), nn.Linear(d // 2, 1))
            for _ in range(n_j_heads)
        ])

        # U: one (real, imag) MLP pair per rep.
        self.u_real_mlps = nn.ModuleList([
            nn.Sequential(nn.Linear(d, d // 2), nn.GELU(), nn.Linear(d // 2, 1))
            for _ in range(n_reps)
        ])
        self.u_imag_mlps = nn.ModuleList([
            nn.Sequential(nn.Linear(d, d // 2), nn.GELU(), nn.Linear(d // 2, 1))
            for _ in range(n_reps)
        ])

    def forward(
        self, z: Tensor, norb: int
    ) -> dict[str, Tensor]:
        B = z.shape[0]
        z_active = z[:, :norb, :norb]
        z_T = z_active.transpose(1, 2)

        n_reps = len(self.u_real_mlps)
        N = z.shape[1]

        j_out = z.new_zeros(B, n_reps, 2, N, N)
        u_real_out = z.new_zeros(B, n_reps, N, N)
        u_imag_out = z.new_zeros(B, n_reps, N, N)

        for rep in range(n_reps):
            for spin in range(2):
                mlp = self.j_mlps[rep * 2 + spin]
                g_pq = mlp(z_active).squeeze(-1)
                g_qp = mlp(z_T).squeeze(-1)
                j_out[:, rep, spin, :norb, :norb] = g_pq + g_qp

            u_real_out[:, rep, :norb, :norb] = self.u_real_mlps[rep](z_active).squeeze(-1)
            u_imag_out[:, rep, :norb, :norb] = self.u_imag_mlps[rep](z_active).squeeze(-1)

        return {
            "j_pred": j_out,
            "u_real_pred": u_real_out,
            "u_imag_pred": u_imag_out,
        }


# ---------------------------------------------------------------------------
# Full pretraining model
# ---------------------------------------------------------------------------


class PretrainingModel(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        self.config = config
        self.tokenizer = Tokenizer(config)
        self.layers = nn.ModuleList([
            TransformerLayer(config) for _ in range(config.num_layers)
        ])
        self.final_norm = nn.LayerNorm(config.embed_dim)
        self.decode_heads = DecodeHeads(config)

        self.apply(self._init_weights)

    def _init_weights(self, module: nn.Module) -> None:
        if isinstance(module, nn.Linear):
            nn.init.trunc_normal_(module.weight, std=0.02)
            if module.bias is not None:
                nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Embedding):
            nn.init.normal_(module.weight, std=0.02)

    def forward(
        self,
        batch: dict[str, Tensor | list],
    ) -> dict[str, Tensor]:
        t2 = batch["t2"]
        B = t2.shape[0]
        nocc = t2.shape[1]
        nvirt = t2.shape[3]
        n_reps = batch["n_reps"][0] if isinstance(batch["n_reps"], list) else batch["n_reps"]
        norb_from_t2 = nocc + nvirt

        norb = norb_from_t2

        orb_energies = batch.get("orb_energies")
        if orb_energies is not None and orb_energies.shape[1] < norb:
            orb_energies = F.pad(orb_energies, (0, norb - orb_energies.shape[1]))

        mo_coeffs = batch.get("mo_coeffs")
        ao_Z = batch.get("ao_Z")
        ao_l = batch.get("ao_l")
        if mo_coeffs is not None and mo_coeffs.shape[1] < norb:
            mo_coeffs = F.pad(mo_coeffs, (0, 0, 0, norb - mo_coeffs.shape[1]))

        x = self.tokenizer(
            t2, nocc, nvirt, n_reps, None,
            orb_energies=orb_energies,
            mo_coeffs=mo_coeffs, ao_Z=ao_Z, ao_l=ao_l,
        )
        t2_bias = self._compute_t2_bias(t2, nocc, nvirt, norb)

        norbs = batch.get("norbs")
        if norbs is not None and isinstance(norbs, Tensor):
            mask = self._compute_batch_mask(norbs, norb, x.device)
        else:
            mask = None

        for layer in self.layers:
            x = layer(x, mask=mask, t2_bias=t2_bias)

        x = self.final_norm(x)
        return self.decode_heads(x, norb)

    def _compute_t2_bias(self, t2: Tensor, nocc: int, nvirt: int, N: int) -> Tensor:
        """Build attention bias from t2 amplitudes.

        For canonical positions (x ∈ occ, a ∈ occ, y ∈ virt):
            bias[x, a, y] = mean_b t2[x, a, y-nocc, b]
        All other positions are zero.
        """
        B = t2.shape[0]
        norb = nocc + nvirt
        bias = t2.new_zeros(B, N, N, N, 1)
        if nocc > 0 and nvirt > 0:
            bias[:, :nocc, :nocc, nocc:norb, 0] = t2.mean(-1)
        return bias

    def _compute_batch_mask(
        self, norbs: Tensor, N: int, device: torch.device
    ) -> Tensor:
        """Per-sample mask: True where any index in the triplet (x,a,y) is padded."""
        idx = torch.arange(N, device=device)
        valid = idx.unsqueeze(0) < norbs.unsqueeze(1).to(device)  # (B, N)
        mask = ~(valid[:, :, None, None] & valid[:, None, :, None] & valid[:, None, None, :])
        return mask  # (B, N, N, N)
