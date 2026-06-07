# Pretrain: Edge Transformer for LUCJ Parameter Prediction

Predicts LUCJ (Locally-Unitary Coupled-Cluster Jastrow) circuit parameters
from CCSD amplitudes using an Edge Transformer architecture adapted from
the (2,1)-GT ("Towards Principled Graph Transformers", NeurIPS 2024).

## Architecture

### Problem Setup

Given a molecular system with $n_\text{orb}$ orbitals ($n_\text{occ}$
occupied + $n_\text{virt}$ virtual), the model takes CCSD amplitudes
$\mathbf{t}_1, \mathbf{t}_2$ as input and predicts the LUCJ diagonal
Coulomb matrices $\mathbf{J}^{(\mu)}$ and orbital rotations
$\mathbf{U}^{(\mu)}$ for $\mu = 1, \ldots, L$ circuit layers.

### Pair-Token Representation

The model operates on **pair tokens** $\mathbf{z}_{pq} \in \mathbb{R}^d$
for all orbital pairs $(p, q)$ with $p, q \in \{0, \ldots, n_\text{orb}-1\}$.
This gives a dense tensor of shape $(B, n_\text{orb}, n_\text{orb}, d)$
throughout the network.

### 1. Tokenizer (Input Embedding)

Each pair token $(p, q)$ is initialized as:

$$\mathbf{z}_{pq}^{(0)} = \text{OrbProj}\bigl([\mathbf{e}_p \| \mathbf{e}_q]\bigr) + \text{PairTypeEmb}(\tau_{pq}) + \text{SliceEnc}(\mathbf{s}_{pq}) + \mathbf{g}$$

| Component | Definition | Implementation |
|:----------|:-----------|:---------------|
| $\mathbf{e}_p$ | See orbital encoding modes below | `--use-hf-energies` flag controls the mode |
| $\text{OrbProj}$ | MLP: $\mathbb{R}^d \to \mathbb{R}^d$ | `Linear(d,d) -> GELU -> Linear(d,d)` |
| $\tau_{pq}$ | Pair type $\in$ {occ-occ, occ-virt, virt-occ, virt-virt} | `nn.Embedding(4, d)` |
| $\mathbf{s}_{pq}$ | $t_2$-slice for pair $(p,q)$, encoded via DeepSets | See below |
| $\mathbf{g}$ | Global: $[n_\text{occ}, n_\text{virt}, L]$ | `Linear(3, d)`, broadcast |

**Orbital encoding modes**:

The per-orbital feature $\mathbf{e}_p \in \mathbb{R}^{d/2}$ is constructed as the
concatenation of a type embedding ($d/4$) and a second embedding ($d/4$) that
varies by mode:

| Mode | Second component of $\mathbf{e}_p$ | Flag | Data source |
|:-----|:-----------------------------------|:-----|:------------|
| **Positional** (default) | $\text{PosEmb}(p)$ — learned per index | (none) | None |
| **HF energies** | $\text{EnergyMLP}(\epsilon_p)$ — Fock eigenvalue | `--use-hf-energies` | `mo_ene_for_qis.dat` |
| **MO coefficients** | $\text{AtomPooling}(\mathbf{C}_p)$ — atom-aware set pooling | `--use-mo-coeffs` | `.in.fchk` |
| **Per-MO graph** (pending) | Equivariant GNN per MO (MoLe-style) | — | `.in.fchk` + coordinates |

**Positional mode**: each orbital gets a learned embedding by index (like
standard Transformer positional encoding). No physical content.

**HF energy mode**: $\epsilon_p$ is the Hartree-Fock orbital energy (Fock
eigenvalue) from the SCF calculation. Encoded via MLP:
`Linear(1, d//4) -> GELU -> Linear(d//4, d//4)`. Orbital energies directly
determine the CC amplitude structure through the energy denominators
$(\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b)^{-1}$.

**MO coefficient mode** (`--use-mo-coeffs`): the richest per-orbital encoding.
Each MO $p$ is defined by its LCAO coefficient vector
$\mathbf{C}_p = (c_{p,1}, \ldots, c_{p,n_\text{basis}})$ in the atomic orbital
basis. Each AO $\mu$ has metadata: nuclear charge $Z_\mu$ (what atom it sits on)
and angular momentum $l_\mu$ (s, p, d, ...). The encoding uses **atom-aware
set pooling**:

$$\mathbf{h}_p = \frac{1}{|\mathcal{B}|}\sum_{\mu \in \mathcal{B}} \phi\bigl(c_{p\mu},\; \text{ZEmb}(Z_\mu),\; \text{lEmb}(l_\mu)\bigr)$$

where $\phi$ is an MLP, $\text{ZEmb}$ and $\text{lEmb}$ are learned embeddings
for nuclear charge and angular momentum, and $\mathcal{B}$ is the set of basis
functions. This handles variable basis sizes naturally (STO-3G has 18--51 AOs
across the dataset) and captures *which kinds of atoms and angular momenta
contribute to each orbital*.

**Per-MO graph embedding** (pending, not yet implemented): following the MoLe
architecture (Thiede et al., arXiv:2602.20232, Feb 2026), each MO would be
embedded as its own molecular graph with equivariant GNN features initialized
from the MO coefficients. This requires 3D atomic coordinates and an
SE(3)-equivariant message-passing network (e.g., MACE). This is the most
physically principled approach — it preserves the spatial structure of each
orbital and respects rotation symmetry — but significantly more complex to
implement. See also: CEONET (King et al., PNAS 2025) for Cartesian equivariant
orbital representations, and bipartite Cholesky graph networks
(arXiv:2605.25268) for integral-factorized orbital encodings.

**$t_2$-slice extraction** for pair $(p, q)$:

| Pair type | Slice from $t_2[i,j,a,b]$ | Shape |
|:----------|:--------------------------|:------|
| occ-occ | $t_2[p, q, :, :]$ | $n_\text{virt} \times n_\text{virt}$ |
| occ-virt | $t_2[p, :, q', :]$ where $q' = q - n_\text{occ}$ | $n_\text{occ} \times n_\text{virt}$ |
| virt-occ | $t_2[:, q, p', :]$ where $p' = p - n_\text{occ}$ | $n_\text{occ} \times n_\text{virt}$ |
| virt-virt | $t_2[:, :, p', q']$ | $n_\text{occ} \times n_\text{occ}$ |

Each slice is flattened and encoded via a **DeepSets** encoder:
$\text{SliceEnc}(\mathbf{s}) = \rho\bigl(\text{mean}(\phi(\mathbf{s}))\bigr)$
where $\phi: \mathbb{R} \to \mathbb{R}^{d/2}$ and $\rho: \mathbb{R}^{d/2} \to \mathbb{R}^d$ are linear layers.

### 2. Edge Attention with $t_2$ Bias

The core attention mechanism uses **factorized pair-to-pair attention**
from the Edge Transformer, with an additive structural bias from $t_2$:

$$A_{(x,a) \to (a,y)} = \frac{\langle \mathbf{z}_{xa} W_Q,\; \mathbf{z}_{ay} W_K \rangle}{\sqrt{d_k}} + \lambda \cdot B_{xay}$$

where:
- The factorized pattern: pair $(x, a)$ attends to pair $(a, y)$ through the
  shared index $a$ (triangular attention).
- $\lambda$ is a **learnable scalar** (initialized to 1.0).
- $B_{xay}$ is a data-dependent bias: for $x, a \in \text{occ}$ and
  $y \in \text{virt}$, $B_{xay} = \text{mean}_b\, t_2[x, a, y', b]$;
  zero otherwise.

The attention computation follows the Edge Transformer einsum pattern:

```
scores = einsum("bxahd, bayhd -> bxayh", Q, K) / sqrt(d_k) + lambda * bias
attn   = softmax(scores, dim=2)  # over shared index a
val    = einsum("bxahd, bayhd -> bxayhd", V_L, V_R)
output = einsum("bxayh, bxayhd -> bxyhd", attn, val)
```

### 3. Transformer Body

Pre-norm Transformer with $L_\text{layers}$ layers:

```
for each layer:
    z' = LayerNorm(z)
    z  = z + EdgeAttention(z', bias=t2_bias)
    z' = LayerNorm(z)
    z  = z + FFN(z')
```

FFN: `Linear(d, 4d) -> GELU -> Dropout -> Linear(4d, d) -> Dropout`

### 4. Decode Heads

From the final pair embeddings $\mathbf{z}_{pq}^{(L)}$, two heads
predict per LUCJ layer $\mu$:

**J head** (symmetric, per spin type):

$$J_\mu^{(\sigma)}[p, q] = g_\sigma^\mu(\mathbf{z}_{pq}) + g_\sigma^\mu(\mathbf{z}_{qp})$$

where $g: \mathbb{R}^d \to \mathbb{R}$ is `Linear(d, d//2) -> GELU -> Linear(d//2, 1)`.
Symmetry is enforced by construction. There are 4 heads total
(2 reps $\times$ 2 spin types: $\alpha\alpha$ and $\alpha\beta$).

**Kappa head** (anti-Hermitian generator of U):

Rather than predicting $U$ directly, the model predicts
$\kappa = \log(U)$, the matrix logarithm. Since $U$ is unitary,
$\kappa$ is anti-Hermitian: its real part is anti-symmetric and its
imaginary part is symmetric. These constraints are enforced by
construction:

$$\kappa_\text{re}^\mu[p,q] = f^\mu(\mathbf{z}_{pq}) - f^\mu(\mathbf{z}_{qp}) \quad \text{(anti-symmetric)}$$
$$\kappa_\text{im}^\mu[p,q] = h^\mu(\mathbf{z}_{pq}) + h^\mu(\mathbf{z}_{qp}) \quad \text{(symmetric)}$$

At inference, $U = \exp(\kappa_\text{re} + i\,\kappa_\text{im})$ is
reconstructed via matrix exponential.

**Why kappa instead of U?** Diagnostic analysis revealed that predicting
$U$ directly caused the model to collapse to near-zero predictions. The
target $U$ matrices are dense complex unitary matrices with entries
spread across $[-0.8, 0.8]$ — predicting all $n_\text{orb}^2$ entries
accurately is extremely hard. The MSE of predicting zeros was ~102 per
sample, and the trained model could not improve beyond this baseline
(U loss stuck at ~131, *worse* than zero). Meanwhile, $\kappa$ provides
structured targets with built-in anti-Hermitian constraints that
regularize the output space and halve the effective number of free
parameters.

See `diagnose_loss.png` for the visualization that motivated this change.

### 5. Loss Function

$$\mathcal{L} = \frac{1}{K} \sum_{k=1}^K \left[ \sum_\mu \bigl\| \hat{J}_\mu^{(k)} - J_\mu^{*(k)} \bigr\|_F^2 \;+\; \alpha \sum_\mu \bigl( \bigl\| \hat{\kappa}_{\text{re},\mu}^{(k)} - \kappa_{\text{re},\mu}^{*(k)} \bigr\|_F^2 + \bigl\| \hat{\kappa}_{\text{im},\mu}^{(k)} - \kappa_{\text{im},\mu}^{*(k)} \bigr\|_F^2 \bigr) \right]$$

- $\alpha$: balancing weight (default 1.0)
- J loss is computed on non-padded entries only
- Kappa loss is computed on non-padded entries only
- U loss is reported for monitoring (via $U = \exp(\hat\kappa)$) but not
  used for gradient computation

## Data Loading

`CCSDAmplitudeDataset` scans the `jobs/` directory for all subdirectories
with both `.done` and `.lucj_done` markers. Each job provides:

**Inputs:**
- `ccsd_t1.dat` $\to$ `t1`: $(n_\text{occ}, n_\text{virt})$
- `ccsd_t2.dat` $\to$ `t2`: $(n_\text{occ}, n_\text{occ}, n_\text{virt}, n_\text{virt})$

**Targets:**
- `lucj_diag_coulomb_mats.npy` $\to$ `j_target`: $(L, 2, n_\text{orb}, n_\text{orb})$
- `kappa_real.npy` $\to$ `kappa_real_target`: $(L, n_\text{orb}, n_\text{orb})$ (precomputed via `precompute_kappa.py`)
- `kappa_imag.npy` $\to$ `kappa_imag_target`: $(L, n_\text{orb}, n_\text{orb})$

**Variable-size handling:** Molecules have $n_\text{orb}$ from 30 to 88.
The `collate_fn` pads all tensors to the maximum size in the batch and
produces boolean masks. A `BucketBatchSampler` groups molecules by similar
$n_\text{orb}$ to minimize padding waste.

## Files

| File | Description |
|:-----|:------------|
| `dataset.py` | `CCSDAmplitudeDataset` with lazy loading, padding collate, norb caching |
| `model.py` | `PretrainingModel`: Tokenizer + EdgeTransformer + DecodeHeads (~1M params at default config) |
| `train.py` | Training loop: AdamW, cosine warmup LR, gradient clipping, checkpointing |
| `plot_curves.py` | Live-updating train/val loss plots from JSON log |

## Usage

### Training

From the `ccsd_amplitudes/` directory:

```bash
# Full training on GPU (run on a machine with GPU, e.g. scai5)
python -m pretrain.train \
    --jobs-dir jobs \
    --epochs 200 \
    --batch-size 4 \
    --lr 1e-4 \
    --embed-dim 128 \
    --num-layers 4 \
    --no-amp \
    --dropout 0.0

# Or via the convenience script (dispatches to GPU server)
bash run_train.sh
```

Key arguments:

| Argument | Default | Description |
|:---------|:--------|:------------|
| `--jobs-dir` | (required) | Path to job directories |
| `--epochs` | 200 | Training epochs |
| `--batch-size` | 32 | Batch size (reduce for large norb) |
| `--lr` | 3e-4 | Peak learning rate |
| `--warmup-steps` | 1000 | Linear LR warmup steps |
| `--alpha` | 1.0 | Weight for U loss term |
| `--embed-dim` | 192 | Model embedding dimension |
| `--num-layers` | 6 | Transformer layers |
| `--num-heads` | 8 | Attention heads |
| `--no-amp` | False | Disable mixed precision (recommended for stability) |
| `--dropout` | 0.1 | Dropout rate (set to 0.0 if NaN issues) |
| `--resume` | None | Path to checkpoint for resuming |

### Monitoring

Training writes a JSON-lines log to `checkpoints/train_log.jsonl`.
Plot train/val curves:

```bash
# One-shot plot
python -m pretrain.plot_curves --log-dir checkpoints --interval 0

# Continuous updates every 30s
python -m pretrain.plot_curves --log-dir checkpoints --interval 30
```

Output: `train_val_curves.png` with three panels (total loss, J loss, U loss).

### Checkpoints

Saved to `checkpoints/`:
- `best.pt`: Best validation loss model
- `last.pt`: Most recent checkpoint
- `epoch_NNNN.pt`: Periodic saves every 10 epochs

Each checkpoint contains: model state_dict, optimizer state, scheduler state,
epoch, global step, and best validation loss. Resume with `--resume checkpoints/last.pt`.

### Inference

Load a trained model and predict LUCJ parameters for a new molecule:

```python
import torch
from pretrain.model import PretrainingModel, ModelConfig
from pretrain.dataset import CCSDAmplitudeDataset

# Load model
config = ModelConfig(embed_dim=128, num_layers=4)
model = PretrainingModel(config)
ckpt = torch.load("checkpoints/best.pt", map_location="cpu")
model.load_state_dict(ckpt["model_state_dict"])
model.eval()

# Load a single sample
ds = CCSDAmplitudeDataset("jobs")
sample = ds[0]
batch = CCSDAmplitudeDataset.collate_fn([sample])

# Predict
with torch.no_grad():
    outputs = model(batch)

j_pred = outputs["j_pred"]       # (1, n_reps, 2, norb, norb)
u_real = outputs["u_real_pred"]  # (1, n_reps, norb, norb)
u_imag = outputs["u_imag_pred"]  # (1, n_reps, norb, norb)

# Reconstruct complex U
U = u_real + 1j * u_imag  # (1, n_reps, norb, norb)
```

## Model Configuration

Default config (`ModelConfig`):

| Parameter | Value | Notes |
|:----------|:------|:------|
| `embed_dim` | 192 | Token embedding dimension $d$ |
| `num_layers` | 6 | Transformer layers |
| `num_heads` | 8 | Attention heads ($d_k = d / H$) |
| `ffn_multiplier` | 4 | FFN hidden dim = $4d$ |
| `dropout` | 0.1 | Dropout rate |
| `attention_dropout` | 0.1 | Attention dropout |
| `max_norb` | 96 | Max orbital count for positional embeddings |
| `n_reps` | 2 | Number of LUCJ circuit layers |

For the current dataset (norb up to 88, batch_size=4), a reduced config
(embed_dim=128, num_layers=4) fits in ~40GB GPU memory.

## References

- Edge Transformer: Mueller et al., "Towards Principled Graph Transformers",
  NeurIPS 2024. arXiv:2401.10119
- Compressed Double Factorization: Lin et al. (2025), arXiv:2511.22476
- ffsim: `UCJOpSpinBalanced.from_t_amplitudes` for LUCJ target generation
