"""Task 2 (Oct 2026): per-molecule direct optimization of LUCJ parameters -- shared pieces.

* Parameter vector x = ffsim UCJOpSpinBalanced parameters (n_reps = 2, square interaction pairs), the
  parameterization of Lin et al.'s NOMAD optimization (src/lucj/quimb_task/lucj_sqd_quimb_task_nomad.py:
  UCJOpSpinBalanced.from_parameters).  Unlike her code, the final orbital rotation is NOT in x: it stays the
  t1 rotation, as for every candidate of this project (the network does not predict it).
  Per layer: norb^2 orbital-rotation parameters (real/imag parts of log U), then J_aa on (p, p+1), then J_ab on
  (p, p).  norb 15 -> 508 parameters, norb 16 -> 574.
* Start points: the optimize=True label (rhf_targets_compressed_small/square_reg0.005) and the network+RL policy
  rl4n29f (runs_ot/energy_tasks/smallval_rl4n29.pkl, made by pretrain.rl.policy_dump from
  rl_runs/grpo_n29_tn/policy_last.pt).
* Exact sampling of the LUCJ state from the GPU CI matrix (pretrain.rl.gpu_energy.LUCJEnergyGPU.state) and
  QSCI / SQD with qiskit-addon-sqd at Lin et al.'s settings (pyscf fixed-space solver instead of Dice).
"""
from __future__ import annotations

import hashlib
import json
import pickle
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HAM_DIR = ROOT / "rhf_hamiltonians"
LABEL_DIR = ROOT / "rhf_targets_compressed_small" / "square_reg0.005"
RL_TASKS = ROOT / "runs_ot" / "energy_tasks" / "smallval_rl4n29.pkl"
RESULTS = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "task2_direct_opt"
LOGS = ROOT / "rl_runs" / "followups" / "task2_direct_opt"
N_REPS = 2

# Lin et al. (arXiv:2511.22476, Sec. "Computational details"; scripts/quimb/*/lucj_compressed_t2_nomad_r24.py)
QSCI_OPT = dict(shots=10_000, samples_per_batch=4000, n_batches=10, max_dim=4000, max_iterations=1,
                symmetrize_spin=True, energy_tol=1e-5, occupancies_tol=1e-3, carryover_threshold=1e-3)
# scoring protocol = Task 1's "lin" protocol (pretrain/followups/task1_sqd_lin.py, seed 0), copied here so that this
# task does not depend on that module while it evolves: her state-vector QSCI (scripts/sqd/*/random_sqd.py; paper:
# 1e5 samples, 10 batches of 4000): ffsim.sample_state_vector(1e6 shots, default_rng(seed)), a uniformly random
# subset of 1e5 (default_rng(12345 + seed)), first iteration of diagonalize_fermionic_hamiltonian, solve_sci with
# spin_sq = 0 (her *_sci variant)
T1_PROTOCOL = dict(shots_raw=1_000_000, shots=100_000, samples_per_batch=4000, n_batches=10, max_dim=4000,
                   max_iterations=1, symmetrize_spin=True, energy_tol=1e-5, occupancies_tol=1e-3,
                   carryover_threshold=1e-3, spin_sq=0.0, seed=0)
NOMAD_PARAMS = ["BB_OUTPUT_TYPE OBJ", "MAX_BB_EVAL 500", "DISPLAY_DEGREE 2", "DISPLAY_ALL_EVAL false",
                "PSD_MADS_OPTIMIZATION True", "PSD_MADS_NB_VAR_IN_SUBPROBLEM 20",
                "PSD_MADS_SUBPROBLEM_MAX_BB_EVAL 20", "PSD_MADS_NB_SUBPROBLEM 4"]


# ----------------------------------------------------------------------------------------- molecules

def load_ham(name: str) -> dict:
    d = np.load(HAM_DIR / f"{name}.npz")
    return dict(one_body=d["one_body"], two_body=d["two_body"], constant=float(d["constant"]),
                norb=int(d["norb"]), nelec=(int(d["nelec_a"]), int(d["nelec_b"])), e_hf=float(d["e_hf"]),
                e_ccsd=float(d["e_ccsd"]))


def load_t1(name: str) -> np.ndarray:
    return np.load(ROOT / "rhf_dataset" / f"{name}.npz")["t1"].astype(np.float64)


def start_uz(name: str, start: str):
    """(U (R,n,n) complex, Z (R,n,n) real) of a start point: 'label' or 'rl4n29f' (or 'rl4n29')."""
    if start == "label":
        lab = np.load(LABEL_DIR / f"{name}.npz")
        return lab["U_re"] + 1j * lab["U_im"], lab["Z"].astype(np.float64)
    tasks = pickle.load(open(RL_TASKS, "rb"))["tasks"]
    for key, n, U, Z, _t1 in tasks:
        if n == name and key[1] == start:
            return np.asarray(U, dtype=complex), np.asarray(Z, dtype=np.float64)
    raise KeyError((name, start))


def corr_pct(E: float, e_hf: float, e_ccsd: float) -> float:
    return 100.0 * (e_hf - E) / (e_hf - e_ccsd)


# ----------------------------------------------------------------------------------------- parameters

def square_pairs(norb: int):
    return [(p, p + 1) for p in range(norb - 1)], [(p, p) for p in range(norb)]


def polar(U: np.ndarray) -> np.ndarray:
    W, _, Vh = np.linalg.svd(U)
    return W @ Vh


def uz_to_x(U: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """(U, Z) in the project's convention -> ffsim parameter vector (no final orbital rotation)."""
    from pretrain.rl.energy import make_ucj_op
    U = polar(np.asarray(U, dtype=complex))
    op = make_ucj_op(np.asarray(Z, dtype=np.float64), U, "square", t1=None)
    return op.to_parameters(interaction_pairs=square_pairs(U.shape[-1]))


def x_to_uz(x: np.ndarray, norb: int, n_reps: int = N_REPS):
    import ffsim
    from pretrain.rl.energy import z_from_op
    op = ffsim.UCJOpSpinBalanced.from_parameters(
        np.asarray(x, dtype=np.float64), norb=norb, n_reps=n_reps, interaction_pairs=square_pairs(norb),
        with_final_orbital_rotation=False)
    return op.orbital_rotations, z_from_op(op, "square")


def n_params(norb: int, n_reps: int = N_REPS) -> int:
    return n_reps * (norb ** 2 + (norb - 1) + norb)


# ----------------------------------------------------------------------------------------- sampling

def sample_ci_matrix(psi, norb: int, k: int, shots: int, seed: int) -> np.ndarray:
    """`shots` configurations from |psi[a, b]|^2 (torch CI matrix on the GPU, alpha index a, beta index b,
    pyscf/ffsim string order), exactly: rows from the marginal, then b | a by an inverse CDF on the sampled rows.
    Returns int64 bitstrings (beta << norb) | alpha  (ffsim BitstringType.INT / qiskit-addon-sqd convention)."""
    import torch
    from pyscf.fci import cistring
    dev = psi.device
    g = torch.Generator(device=dev)
    g.manual_seed(int(seed))
    dim = psi.shape[0]
    # row marginal in float64, in row blocks (no full-size temporary)
    pr = torch.empty(dim, device=dev, dtype=torch.float64)
    B = max(1, int(2e8 // (16 * dim)))
    for r0 in range(0, dim, B):
        blk = psi[r0:r0 + B]
        pr[r0:r0 + B] = (blk.real.double().square() + blk.imag.double().square()).sum(1)
    tot = float(pr.sum())
    rows = torch.multinomial(pr / tot, shots, replacement=True, generator=g)
    u = torch.rand(shots, device=dev, dtype=torch.float64, generator=g)
    urow, inv = torch.unique(rows, return_inverse=True)
    cols = torch.empty(shots, device=dev, dtype=torch.int64)
    CH = max(1, int(1e8 // (16 * dim)))                  # unique rows per chunk
    for c0 in range(0, len(urow), CH):
        sel = urow[c0:c0 + CH]
        blk = psi[sel]
        w = blk.real.double().square() + blk.imag.double().square()
        cdf = torch.cumsum(w, 1)
        cdf /= cdf[:, -1:].clone()
        m = (inv >= c0) & (inv < c0 + len(sel))
        idx = torch.nonzero(m).squeeze(1)
        rr = inv[idx] - c0
        # flattened search: row r occupies (r, r + 1]; cdf + r is globally sorted
        flat = (cdf + torch.arange(len(sel), device=dev, dtype=torch.float64)[:, None]).reshape(-1)
        pos = torch.searchsorted(flat, rr.double() + u[idx])
        cols[idx] = torch.clamp(pos - rr * dim, 0, dim - 1)
        del blk, w, cdf, flat
    strs = np.asarray(cistring.make_strings(range(norb), k), dtype=np.int64)
    a = rows.cpu().numpy()
    b = cols.cpu().numpy()
    return strs[a] | (strs[b] << norb)


# ----------------------------------------------------------------------------------------- QSCI / SQD

def task1_ci_strings(samples_int, norb: int, nelec, rng, seed: int = 0):
    """CI strings of the 10 batches of the Task-1 "lin" protocol from 1e6 raw samples (identical to
    task1_sqd_lin.lin_ci_strings with PROTOCOLS["lin"] and subset_rng = default_rng(12345 + seed)).
    `rng` is the generator that drew the samples (its state continues into the batch subsampling)."""
    from qiskit.primitives import BitArray
    from qiskit_addon_sqd.counts import bit_array_to_arrays
    from qiskit_addon_sqd.fermion import _LoopConfig, _prepare_ci_strings
    P = T1_PROTOCOL
    s = np.asarray(samples_int, dtype=np.int64)
    array = ((s[:, None] >> np.arange(2 * norb - 1, -1, -1)[None, :]) & 1).astype(bool)
    chk = BitArray.from_samples(list(map(int, s[:1000])), num_bits=2 * norb).to_bool_array()
    assert np.array_equal(chk, array[:1000])
    keep = np.random.default_rng(12345 + seed).choice(np.arange(0, array.shape[0]), size=P["shots"], replace=False)
    bit_array = BitArray.from_bool_array(array[keep])
    raw_bitstrings, raw_probs = bit_array_to_arrays(bit_array)
    empty = np.array([], dtype=np.int64)
    cfg = _LoopConfig(raw_bitstrings=raw_bitstrings, raw_probs=raw_probs, n_alpha=nelec[0], n_beta=nelec[1],
                      samples_per_batch=P["samples_per_batch"], num_batches=P["n_batches"], norb=norb,
                      symmetrize_spin=P["symmetrize_spin"], include_a=np.unique(np.array([], dtype=int)),
                      include_b=np.unique(np.array([], dtype=int)), max_dim_a=P["max_dim"], max_dim_b=P["max_dim"],
                      energy_tol=P["energy_tol"], occupancies_tol=P["occupancies_tol"],
                      carryover_threshold=P["carryover_threshold"], rng=np.random.default_rng(rng))
    ci = _prepare_ci_strings(cfg, None, empty, empty)
    return ci, int(raw_bitstrings.shape[0])


class CachedSolver:
    """qiskit-addon-sqd solve_sci_batch (pyscf kernel_fixed_space) that diagonalizes identical subspaces once
    (with fewer than samples_per_batch distinct valid bitstrings every batch is the full set; results are
    identical, only the repeats are skipped)."""

    def __init__(self):
        self.n_solved = 0
        self.t_solve = 0.0

    def __call__(self, ci_strings, one_body, two_body, norb, nelec):
        from qiskit_addon_sqd.fermion import solve_sci
        out, memo = [], {}
        for sa, sb in ci_strings:
            key = hashlib.sha1(np.ascontiguousarray(sa).tobytes() + b"|" + np.ascontiguousarray(sb).tobytes()).hexdigest()
            if key not in memo:
                t0 = time.time()
                memo[key] = solve_sci((sa, sb), one_body, two_body, norb=norb, nelec=nelec)
                self.t_solve += time.time() - t0
                self.n_solved += 1
            out.append(memo[key])
        return out


def qsci_energy(ham: dict, bitstrings: np.ndarray, settings: dict, rng, keep_best: bool = False) -> dict:
    """QSCI as in Lin et al.: diagonalize_fermionic_hamiltonian(..., max_iterations=1, symmetrize_spin) on the
    given samples.  Returns the per-batch energies (+ constant), subspace dims and the returned (lowest) energy.
    Run it only in a process without torch (task2_sci_worker): see the note there."""
    from qiskit.primitives import BitArray
    from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian
    norb, nelec = ham["norb"], ham["nelec"]
    t0 = time.time()
    bit_array = BitArray.from_samples([int(v) for v in bitstrings], num_bits=2 * norb)
    hist = []

    def cb(results):
        hist.append([(float(r.energy + ham["constant"]), tuple(int(s) for s in r.sci_state.amplitudes.shape))
                     for r in results])

    solver = CachedSolver()
    res = diagonalize_fermionic_hamiltonian(
        ham["one_body"], ham["two_body"], bit_array, samples_per_batch=settings["samples_per_batch"],
        norb=norb, nelec=nelec, num_batches=settings["n_batches"], energy_tol=settings["energy_tol"],
        occupancies_tol=settings["occupancies_tol"], max_iterations=settings["max_iterations"],
        sci_solver=solver, symmetrize_spin=settings["symmetrize_spin"],
        carryover_threshold=settings["carryover_threshold"], seed=rng, max_dim=settings["max_dim"], callback=cb)
    last = hist[-1]
    Es = [e for e, _ in last]
    out = dict(E=float(res.energy + ham["constant"]), E_batches=Es, E_mean=float(np.mean(Es)),
               E_min=float(np.min(Es)), E_max=float(np.max(Es)), dims=[d for _, d in last],
               n_unique=int(len(np.unique(bitstrings))), n_solved=solver.n_solved, t_solve=solver.t_solve,
               t=time.time() - t0)
    if keep_best:
        out["best"] = res
    return out


def dumpj(obj, path: Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    json.dump(obj, open(tmp, "w"), indent=1)
    tmp.replace(path)
