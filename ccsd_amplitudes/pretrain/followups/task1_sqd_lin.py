#!/usr/bin/env python3
"""Task 1 (Oct 2026): noiseless QSCI ("SQD" without configuration recovery) of LUCJ candidates with the protocol of
Lin et al., arXiv:2511.22476.  Question: does a better LUCJ variational energy give a better QSCI energy?

Primary protocol "lin" (paper, "State vector simulation and QSCI"; her scripts/sqd/fe2s2_30e20o/lucj_*_t2.py and
n2_cc-pvdz_10e26o/random_sqd.py: shots=100_000, n_batches=10, samples_per_batch=max_dim=4000; task code
src/lucj/sqd_energy_task/lucj_compressed_t2_task{,_sci}.py):
  1. exact LUCJ state; ffsim.sample_state_vector(shots=1_000_000, seed=rng), rng = default_rng(seed) (her entropy=0
     is our seed 0);
  2. keep a uniformly random subset of 100_000 of these samples (she uses the unseeded global np.random; we seed it);
  3. qiskit_addon_sqd.fermion.diagonalize_fermionic_hamiltonian(..., samples_per_batch=4000, num_batches=10,
     max_dim=4000, max_iterations=1, symmetrize_spin=True, energy_tol=1e-5, occupancies_tol=1e-3,
     carryover_threshold=1e-3, seed=rng): one iteration, i.e. no configuration recovery;
  4. solver: PySCF selected CI, qiskit_addon_sqd.fermion.solve_sci_batch with spin_sq=0.0 (her *_sci variant;
     her Dice variant needs the Dice binary, which is not installed here);
  5. report the mean / min / max of the 10 batch energies (her figures); her code's returned value is the min.
Secondary protocol "n2631g": her script for the only system of our size, scripts/sqd/n2_6-31g_10e16o/lucj_*.py:
  shots=1_000_000, one batch with every sampled configuration, max_dim=None (same solver, one iteration).
The SQD loop is split into stages so the diagonalizations can run anywhere.  `sample` builds the CI strings with
the library's own _prepare_ci_strings (the first iteration of diagonalize_fermionic_hamiltonian; that code is
identical in qiskit-addon-sqd 0.12.0 (her lock file) and 0.13.1 (ours)); `selftest` checks that the staged
pipeline reproduces her one-call diagonalize_fermionic_hamiltonian bit for bit.  Identical subspaces (all batches
are identical when 10^5 samples hold fewer than 4,000 distinct configurations) are diagonalized once.

Stages:
  candidates  (CPU)  policy_dump (label, pre4, rl4L, rl4n29f) + truncated CCSD init [+ Lin-style compressed DF]
  sample      (GPU)  exact state (gpu_energy.LUCJEnergyGPU), exact LUCJ energy, samples, CI strings per protocol
  diag        (CPU)  PySCF selected-CI diagonalizations; one claim per sample file, several workers can share
  fci         (CPU)  PySCF FCI references (direct_spin0)
  report             summary JSON (tables: task1_sqd_report.py)
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import pickle
import socket
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

TASK = "task1_sqd_vs_lin"
WORK = ROOT / "runs_ot" / "energy_tasks" / TASK          # bulky intermediates (gitignored)
RES = ROOT / "pretrain" / "opt_true" / "results" / "followups" / TASK
LOGS = ROOT / "rl_runs" / "followups" / TASK
HAM_DIR = ROOT / "rhf_hamiltonians"

# Lin et al. settings common to both protocols
SQD = dict(max_iterations=1, symmetrize_spin=True, energy_tol=1e-5, occupancies_tol=1e-3, carryover_threshold=1e-3,
           spin_sq=0.0, shots_raw=1_000_000)
PROTOCOLS = {
    # paper / fe2s2 + cc-pVDZ scripts
    "lin": dict(shots=100_000, samples_per_batch=4000, n_batches=10, max_dim=4000),
    # n2_6-31g_10e16o scripts (16 orbitals, like ours)
    "n2631g": dict(shots=1_000_000, samples_per_batch=1_000_000, n_batches=1, max_dim=None),
}

CAND_ORDER = ["truncated", "cdf_lin", "label", "pre4", "rl4L", "rl4n29f"]
CAND_LABEL = {"truncated": "truncated CCSD (n_reps=2, square)",
              "cdf_lin": "compressed DF, Lin et al. settings (no reg.)",
              "label": "optimize=True label (reg. 0.005, <=500 it.)",
              "pre4": "pretrained network, 4 recycles",
              "rl4L": "network + RL (norb 15-18, step 25)",
              "rl4n29f": "network + RL (+ n29 TN RL, step 30)"}


def set_threads(n: int):
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS",
              "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[k] = str(n)
    os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
    os.environ.setdefault("JAX_PLATFORMS", "cpu")


def names_from(path):
    return [ln.strip() for ln in open(ROOT / path) if ln.strip()]


def load_ham(name):
    d = np.load(HAM_DIR / f"{name}.npz")
    return dict(h1=d["one_body"].astype(np.float64), h2=d["two_body"].astype(np.float64),
                const=float(d["constant"]), norb=int(d["norb"]), nelec=(int(d["nelec_a"]), int(d["nelec_b"])),
                e_hf=float(d["e_hf"]), e_ccsd=float(d["e_ccsd"]))


# ===================================================================== candidates
def truncated_uz(t1, t2, n_reps=2, connectivity="square"):
    """Lin et al.'s 'truncated' init: ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, t1=t1, n_reps, square pairs)
    (her load_operator, truncated branch).  Returns (U, Z) in this project's convention (make_ucj_op)."""
    import ffsim
    from ffsim.variational.util import interaction_pairs_spin_balanced

    from pretrain.rl.energy import z_from_op
    norb = t2.shape[0] + t2.shape[2]
    op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2, t1=t1, n_reps=n_reps, interaction_pairs=interaction_pairs_spin_balanced(connectivity, norb))
    return np.asarray(op.orbital_rotations), z_from_op(op, connectivity), op


def cdf_lin_uz(t1, t2, n_reps=2, connectivity="square", maxiter=100, begin_reps=20, step=2):
    """Lin et al.'s 'compressed' init (paper, Computational details): compressed double factorization by L-BFGS-B,
    maxiter 100, multi-stage from 20 repetitions in steps of 2, no regularization (her QSCI choice), allowed
    Coulomb entries = union of the square pairs.  ffsim's implementation (she contributed it to ffsim)."""
    import ffsim
    from ffsim.variational.util import interaction_pairs_spin_balanced

    from pretrain.rl.energy import z_from_op
    norb = t2.shape[0] + t2.shape[2]
    op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2, t1=t1, n_reps=n_reps, interaction_pairs=interaction_pairs_spin_balanced(connectivity, norb),
        optimize=True, options=dict(maxiter=maxiter), multi_stage_start=begin_reps, multi_stage_step=step)
    return np.asarray(op.orbital_rotations), z_from_op(op, connectivity), op


def cmd_candidates(a):
    set_threads(a.threads)
    WORK.mkdir(parents=True, exist_ok=True)
    out = WORK / f"cands_{a.tag}.pkl"
    pd_out = WORK / f"policy_dump_{a.tag}.pkl"
    if not pd_out.exists():
        cmd = [sys.executable, "-m", "pretrain.rl.policy_dump", "--names-file", a.names_file,
               "--policy", "init4:runs_ot/slotall_T4/best.pt", "--tag", "pre4",
               "--policy", "rl_runs/grpo_large_rl4/policy_best.pt", "--tag", "rl4L",
               "--policy", "rl_runs/grpo_n29_tn/policy_last.pt", "--tag", "rl4n29f",
               "--labels-dir", a.labels_dir, "--device", "cpu", "--out", str(pd_out.relative_to(ROOT))]
        print(" ".join(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
    pd = pickle.load(open(pd_out, "rb"))
    tasks, meta = list(pd["tasks"]), pd["meta"]
    names = names_from(a.names_file)
    for n in names:
        d = np.load(ROOT / "rhf_dataset" / f"{n}.npz")
        t1, t2 = d["t1"].astype(np.float64), d["t2"].astype(np.float64)
        t0 = time.time()
        U, Z, _ = truncated_uz(t1, t2)
        tasks.append(((n, "truncated"), n, U, Z, t1))
        meta[n].setdefault("time", {})["truncated"] = time.time() - t0
        if a.with_cdf:
            t0 = time.time()
            U, Z, _ = cdf_lin_uz(t1, t2)
            tasks.append(((n, "cdf_lin"), n, U, Z, t1))
            meta[n]["time"]["cdf_lin"] = time.time() - t0
            print(f"  {n}: cdf_lin {time.time() - t0:.1f}s", flush=True)
    pickle.dump({"tasks": tasks, "meta": meta, "args": vars(a)}, open(out, "wb"))
    print(f"{len(tasks)} candidate tasks -> {out}")


# ===================================================================== sampling (GPU)
def lin_ci_strings(samples_int, norb, nelec, rng, *, shots, spb, n_batches, max_dim, subset_rng):
    """Her post-sampling steps up to the solver call: BitArray of the raw samples, uniformly random subset of
    `shots` rows, then the first iteration of diagonalize_fermionic_hamiltonian (postselection, subsample,
    symmetrize, max_dim) with qiskit_addon_sqd's own helpers.  Returns (ci_strings, n_unique_in_subset)."""
    from qiskit.primitives import BitArray
    from qiskit_addon_sqd.counts import bit_array_to_arrays
    from qiskit_addon_sqd.fermion import _LoopConfig, _prepare_ci_strings

    s = np.asarray(samples_int, dtype=np.int64)
    # == BitArray.from_samples(samples, num_bits=2 * norb).to_bool_array() (big-endian rows; checked on the first
    # 1000 samples every call), vectorized because from_samples loops over Python ints
    array = ((s[:, None] >> np.arange(2 * norb - 1, -1, -1)[None, :]) & 1).astype(bool)
    chk = BitArray.from_samples(list(map(int, s[:1000])), num_bits=2 * norb).to_bool_array()
    assert np.array_equal(chk, array[:1000])
    keep = subset_rng.choice(np.arange(0, array.shape[0]), size=shots, replace=False)
    bit_array = BitArray.from_bool_array(array[keep])
    raw_bitstrings, raw_probs = bit_array_to_arrays(bit_array)
    empty = np.array([], dtype=np.int64)
    cfg = _LoopConfig(raw_bitstrings=raw_bitstrings, raw_probs=raw_probs, n_alpha=nelec[0], n_beta=nelec[1],
                      samples_per_batch=spb, num_batches=n_batches, norb=norb,
                      symmetrize_spin=SQD["symmetrize_spin"], include_a=np.unique(np.array([], dtype=int)),
                      include_b=np.unique(np.array([], dtype=int)), max_dim_a=max_dim, max_dim_b=max_dim,
                      energy_tol=SQD["energy_tol"], occupancies_tol=SQD["occupancies_tol"],
                      carryover_threshold=SQD["carryover_threshold"], rng=np.random.default_rng(rng))
    ci = _prepare_ci_strings(cfg, None, empty, empty)
    return ci, int(raw_bitstrings.shape[0])


def sample_file(name, cand, seed):
    return WORK / "samples" / f"{name}__{cand}__s{seed}.npz"


def cmd_sample(a):
    set_threads(a.threads)
    import ffsim
    import torch
    torch.set_num_threads(a.threads)
    from pretrain.rl.gpu_energy import LUCJEnergyGPU, corr_frac
    tasks = pickle.load(open(a.cands, "rb"))["tasks"]
    only = set(a.only_names or [])
    by_mol = {}
    for key, name, U, Z, t1 in tasks:
        if only and name not in only:
            continue
        if a.cands_only and key[1] not in a.cands_only:
            continue
        by_mol.setdefault(name, []).append((key[1], U, Z, t1))
    (WORK / "samples").mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    dev = torch.device("cuda", 0)
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    pool = ProcessPoolExecutor(max_workers=a.ci_workers, mp_context=mp.get_context("spawn"))   # no CUDA in children
    for name, items in by_mol.items():
        todo = [it for it in items if any(not sample_file(name, it[0], s).exists() for s in a.seeds)]
        if not todo:
            continue
        h = load_ham(name)
        eng = LUCJEnergyGPU.from_npz(name, ham_dir=HAM_DIR, device=dev, dtype=torch.complex64,
                                     max_mem_gb=a.max_mem_gb)
        for cand, U, Z, t1 in todo:
            t0 = time.time()
            E = eng.energy(U, Z, t1=t1, connectivity="square")
            t_e = time.time() - t0
            psi = eng.psi.detach().to("cpu").numpy().astype(np.complex128).reshape(-1)
            nrm = float(np.linalg.norm(psi))
            psi /= nrm
            p = np.abs(psi) ** 2
            ent = float(-(p[p > 0] * np.log(p[p > 0])).sum())
            p_hf = float(p[0])            # HF = lowest k orbitals occupied = address 0 in ffsim order
            del p
            pend = []
            for seed in a.seeds:
                if sample_file(name, cand, seed).exists():
                    continue
                t0 = time.time()
                rng = np.random.default_rng(seed)
                samples = ffsim.sample_state_vector(psi, norb=h["norb"], nelec=h["nelec"],
                                                    shots=SQD["shots_raw"], seed=rng,
                                                    bitstring_type=ffsim.BitstringType.INT)
                t_s = time.time() - t0
                meta = dict(name=name, cand=cand, seed=seed, e_lucj=E,
                            corr_frac_ccsd=corr_frac(E, h["e_hf"], h["e_ccsd"]), norm=nrm, entropy=ent,
                            p_hf=p_hf, t_energy=t_e, t_sample=t_s, gpu=torch.cuda.get_device_name(dev),
                            host=socket.gethostname())
                pend.append(pool.submit(postprocess_seed, np.asarray(samples, dtype=np.int64), rng, meta,
                                        h["norb"], h["nelec"]))
            for fut in pend:
                rec = fut.result()
                with open(RES / "lucj_states.jsonl", "a") as f:
                    f.write(json.dumps(rec) + "\n")
                lp, np_ = rec["protocols"]["lin"], rec["protocols"]["n2631g"]
                print(f"  {name:20s} {cand:10s} s{rec['seed']} E {E:.8f} corr {100 * rec['corr_frac_ccsd']:6.2f}%  "
                      f"p_HF {p_hf:.4f}  uniq 1e5 {lp['n_unique_subset']:6d} / 1e6 {rec['n_unique_1M']:6d}  "
                      f"lin dims {min(lp['dims'])}-{max(lp['dims'])} ({lp['n_distinct_subspaces']} distinct)  "
                      f"n2631g dim {np_['dims'][0]}  (E {t_e:.1f}s, sample {rec['t_sample']:.1f}s, "
                      f"ci {rec['t_ci']:.1f}s)", flush=True)
            del psi
        eng.release()
        del eng
        torch.cuda.empty_cache()
    pool.shutdown()


def postprocess_seed(samples, rng, meta, norb, nelec):
    """CI strings of every protocol from one 10^6-sample draw (runs in a worker process); saves the sample file."""
    set_threads(1)
    t0 = time.time()
    seed = meta["seed"]
    out, rec_p = {}, {}
    for pname, pr in PROTOCOLS.items():
        ci, nu = lin_ci_strings(samples, norb, nelec, copy.deepcopy(rng), shots=pr["shots"],
                                spb=pr["samples_per_batch"], n_batches=pr["n_batches"], max_dim=pr["max_dim"],
                                subset_rng=np.random.default_rng(12345 + seed))
        assert all(np.array_equal(sa, sb) for sa, sb in ci)
        for b, (sa, _) in enumerate(ci):
            out[f"{pname}_ci_{b}"] = sa
        rec_p[pname] = dict(n_unique_subset=nu, dims=[int(len(sa)) for sa, _ in ci],
                            n_distinct_subspaces=len({hashlib.sha1(sa.tobytes()).hexdigest() for sa, _ in ci}))
    u, cnt = np.unique(samples, return_counts=True)
    dst = sample_file(meta["name"], meta["cand"], seed)
    part = dst.with_name(dst.name + ".part")
    with open(part, "wb") as fh:                         # atomic: diag workers glob *.npz
        np.savez_compressed(fh, uniq_1M=u, counts_1M=cnt, e_lucj=meta["e_lucj"], norm=meta["norm"],
                            entropy=meta["entropy"], p_hf=meta["p_hf"], norb=norb, nelec=np.array(nelec),
                            protocols=json.dumps(rec_p), **out)
    os.replace(part, dst)
    return dict(meta, n_unique_1M=int(len(u)), protocols=rec_p, t_ci=time.time() - t0)


# ===================================================================== diagonalizations (CPU)
def solve_strings(strs, h, threads):
    from pyscf import lib
    from qiskit_addon_sqd.fermion import solve_sci
    lib.num_threads(threads)
    t0 = time.time()
    res = solve_sci((strs, strs), h["h1"], h["h2"], norb=h["norb"], nelec=h["nelec"], spin_sq=SQD["spin_sq"])
    t_solve = time.time() - t0
    t0 = time.time()
    s2 = float(res.sci_state.spin_square())
    t_s2 = time.time() - t0
    amp = res.sci_state.amplitudes
    hf = (1 << h["nelec"][0]) - 1
    return dict(energy=float(res.energy + h["const"]), dim_a=int(len(strs)), subspace_dim=int(len(strs)) ** 2,
                spin_square=s2, c_hf=float(abs(amp[0, 0])) if int(strs[0]) == hf else None,
                t_solve=t_solve, t_s2=t_s2)


def topk_strings(z, norb, nelec, k):
    """Protocol "top<k>": her n2631g protocol (every configuration of the 10^6 samples, one batch) with max_dim=k,
    i.e. the k symmetrized single-spin strings that the library ranks first (by the number of distinct sampled
    configurations containing them).  Rebuilt from the stored unique samples; with k=None it reproduces n2631g."""
    from qiskit_addon_sqd.fermion import _LoopConfig, _prepare_ci_strings
    u, cnt = z["uniq_1M"].astype(np.int64), z["counts_1M"]
    bits = ((u[:, None] >> np.arange(2 * norb - 1, -1, -1)[None, :]) & 1).astype(bool)
    empty = np.array([], dtype=np.int64)
    cfg = _LoopConfig(raw_bitstrings=bits, raw_probs=cnt / cnt.sum(), n_alpha=nelec[0], n_beta=nelec[1],
                      samples_per_batch=SQD["shots_raw"], num_batches=1, norb=norb, symmetrize_spin=True,
                      include_a=np.unique(np.array([], dtype=int)), include_b=np.unique(np.array([], dtype=int)),
                      max_dim_a=k, max_dim_b=k, energy_tol=SQD["energy_tol"], occupancies_tol=SQD["occupancies_tol"],
                      carryover_threshold=SQD["carryover_threshold"], rng=np.random.default_rng(0))
    return _prepare_ci_strings(cfg, None, empty, empty)[0][0]


def cmd_diag(a):
    """Diagonalize every sample file (all protocols); with --idle-exit-s, keep polling for new sample files."""
    set_threads(a.threads)
    import signal
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))     # run the finally blocks (claim release) on kill
    t_idle = time.time()
    while True:
        n = diag_pass(a)
        if n:
            t_idle = time.time()
        if a.stop_file and Path(a.stop_file).exists():
            break
        if not a.idle_exit_s or time.time() - t_idle > a.idle_exit_s:
            break
        time.sleep(30)


def diag_out(stem, proto):
    return WORK / "diag" / f"{stem}__{proto}.json"


def n_batches_of(proto):
    return PROTOCOLS[proto]["n_batches"] if proto in PROTOCOLS else 1


def diag_pass(a):
    """One sweep over (sample file, protocol) units; a unit is claimed by one worker (O_EXCL claim file) and its
    batch results are saved after every batch, so an interrupted unit resumes where it stopped."""
    if getattr(a, "names_file", None):
        a.only_names = sorted(set(a.only_names or []) | set(names_from(a.names_file)))
    ddir = WORK / "diag"
    cdir = ddir / "claims"
    cdir.mkdir(parents=True, exist_ok=True)
    protos = a.protocols or list(PROTOCOLS)
    files = sorted((WORK / "samples").glob("*.npz"))
    units = []
    for f in files:
        name, cand, s = f.stem.split("__")
        seed = int(s[1:])
        if time.time() - f.stat().st_mtime < 60:         # may still be being written by the sampler
            continue
        if a.only_names and name not in a.only_names:
            continue
        if a.only_cands and cand not in a.only_cands:
            continue
        if a.skip_cands and cand in a.skip_cands:
            continue
        if a.seeds and seed not in a.seeds:
            continue
        for pi, p in enumerate(protos):
            units.append(((pi, seed, name, CAND_ORDER.index(cand) if cand in CAND_ORDER else 99), f, p))
    units.sort(key=lambda u: u[0])
    t_start, n_solved = time.time(), 0
    for _, f, p in units:
        out = diag_out(f.stem, p)
        nb = n_batches_of(p)
        done = json.load(open(out)) if out.exists() else {"batches": []}
        if len(done["batches"]) >= nb:
            continue
        claim = cdir / f"{f.stem}__{p}.claim"
        try:
            fd = os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        with os.fdopen(fd, "w") as fh:
            fh.write(json.dumps(dict(host=socket.gethostname(), pid=os.getpid(), t=time.time())))
        try:
            done = json.load(open(out)) if out.exists() else {"batches": []}
            name = f.stem.split("__")[0]
            h = load_ham(name)
            z = np.load(f)
            cache = {r["sha1"]: r for r in done["batches"]}
            for b in range(len(done["batches"]), nb):
                strs = z[f"{p}_ci_{b}"] if p in PROTOCOLS else topk_strings(z, h["norb"], h["nelec"], int(p[3:]))
                sha = hashlib.sha1(strs.tobytes()).hexdigest()
                if sha in cache:
                    r = dict(cache[sha], batch=b, dedup=True)
                else:
                    r = dict(solve_strings(strs, h, a.threads), batch=b, sha1=sha, dedup=False,
                             threads=a.threads, host=socket.gethostname(), finished=time.time())
                    cache[sha] = r
                    n_solved += 1
                    print(f"  {f.stem:36s} {p:6s} b{b}  E {r['energy']:.8f}  dim {r['dim_a']}^2  "
                          f"S2 {r['spin_square']:.1e}  {r['t_solve']:.1f}s  [{n_solved} solved, "
                          f"{time.time() - t_start:.0f}s]", flush=True)
                done["batches"].append(r)
                tmp = out.with_name(out.name + ".tmp")
                json.dump(done, open(tmp, "w"))
                os.replace(tmp, out)
        finally:
            claim.unlink(missing_ok=True)
        if a.stop_file and Path(a.stop_file).exists():
            print("stop file present, exiting", flush=True)
            break
    return n_solved


# ===================================================================== FCI references (CPU)
def cmd_fci(a):
    set_threads(a.threads)
    from pyscf import fci, lib
    lib.num_threads(a.threads)
    RES.mkdir(parents=True, exist_ok=True)
    path = RES / "fci_refs.json"
    for name in (a.names or names_from(a.names_file)):
        refs = json.load(open(path)) if path.exists() else {}
        if name in refs and not a.redo:
            continue
        h = load_ham(name)
        solver = fci.direct_spin0.FCI()
        solver.max_memory = a.max_memory_mb
        solver.conv_tol = a.conv_tol
        solver.max_cycle = 200
        solver.verbose = 5 if a.verbose else 0
        t0 = time.time()
        e, civec = solver.kernel(h["h1"], h["h2"], h["norb"], h["nelec"], ecore=h["const"])
        dt = time.time() - t0
        ss, _ = solver.spin_square(civec, h["norb"], h["nelec"])
        refs = json.load(open(path)) if path.exists() else {}
        refs[name] = dict(e_fci=float(e), s2=float(ss), conv_tol=a.conv_tol, t_s=dt, threads=a.threads,
                          solver="pyscf.fci.direct_spin0", converged=bool(np.all(solver.converged)),
                          dim=int(civec.size), c_hf=float(abs(civec.reshape(-1)[0])), host=socket.gethostname())
        json.dump(refs, open(path, "w"), indent=1)
        print(f"  {name:20s} E_FCI {e:.10f}  S2 {ss:.1e}  {dt:.0f}s  dim {civec.size}", flush=True)
        del civec


# ===================================================================== self-test of the staged pipeline
def cmd_selftest(a):
    """Small molecule: staged pipeline (lin_ci_strings + solve_sci per batch) == her one-call
    diagonalize_fermionic_hamiltonian with the same generator.  Also ffsim vs GPU state (with --gpu-check)."""
    set_threads(a.threads)
    from functools import partial

    import ffsim
    from pyscf import cc, gto, scf
    from qiskit.primitives import BitArray
    from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci, solve_sci_batch

    from pretrain.rl.energy import make_ucj_op
    mol = gto.M(atom="N 0 0 0; N 0 0 1.2", basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    active = list(range(2, 12))
    md = ffsim.MolecularData.from_scf(mf, active_space=active)
    norb, nelec = md.norb, md.nelec
    mycc = cc.CCSD(mf, frozen=[i for i in range(mol.nao) if i not in active]).run()
    U, Z, op = truncated_uz(mycc.t1, mycc.t2)
    op2 = make_ucj_op(Z, U, "square", t1=mycc.t1)
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, nelec), op, norb=norb, nelec=nelec)
    psi2 = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, nelec), op2, norb=norb, nelec=nelec)
    print(f"truncated op via (U,Z) round trip: max|dpsi| = {np.abs(psi - psi2).max():.2e}")
    ham = md.hamiltonian
    if a.gpu_check:
        import torch

        from pretrain.rl.gpu_energy import LUCJEnergyGPU
        for dt in (torch.complex128, torch.complex64):
            eng = LUCJEnergyGPU(ham.one_body_tensor, ham.two_body_tensor, ham.constant, norb, nelec,
                                device="cuda:0", dtype=dt)
            E = eng.energy(U, Z, t1=mycc.t1)
            pg = eng.psi.cpu().numpy().reshape(-1).astype(np.complex128)
            e_ff = float(np.vdot(psi, ffsim.linear_operator(ham, norb=norb, nelec=nelec) @ psi).real)
            print(f"GPU ({dt}) vs ffsim state: |<ffsim|gpu>| = {abs(np.vdot(psi, pg)):.12f}, "
                  f"max|psi_gpu - psi_ffsim| = {np.abs(pg - psi).max():.2e};  E_gpu {E:.10f}  E_ffsim {e_ff:.10f}")
    # a spread-out state (random LUCJ, as her random_sqd baseline) so that subsampling and max_dim both matter
    rop = ffsim.random.random_ucj_op_spin_balanced(norb, n_reps=2, with_final_orbital_rotation=True, seed=7)
    psi = ffsim.apply_unitary(ffsim.hartree_fock_state(norb, nelec), rop, norb=norb, nelec=nelec)
    spb, nb, md_ = 300, 4, 200
    # (a) her one-call path
    rng = np.random.default_rng(0)
    samples = ffsim.sample_state_vector(psi, norb=norb, nelec=nelec, shots=20000, seed=rng,
                                        bitstring_type=ffsim.BitstringType.INT)
    bit_array = BitArray.from_samples(samples, num_bits=2 * norb)
    array = bit_array.to_bool_array()
    keep = np.random.default_rng(12345).choice(np.arange(0, array.shape[0]), size=5000, replace=False)
    ba = BitArray.from_bool_array(array[keep])
    hist = []
    res = diagonalize_fermionic_hamiltonian(
        ham.one_body_tensor, ham.two_body_tensor, ba, samples_per_batch=spb, norb=norb, nelec=nelec,
        num_batches=nb, energy_tol=SQD["energy_tol"], occupancies_tol=SQD["occupancies_tol"], max_iterations=1,
        sci_solver=partial(solve_sci_batch, spin_sq=0.0), symmetrize_spin=True,
        carryover_threshold=SQD["carryover_threshold"], seed=rng,
        callback=lambda rs: hist.append([(r.energy, r.sci_state.amplitudes.shape) for r in rs]), max_dim=md_)
    # (b) staged path
    rng = np.random.default_rng(0)
    samples_b = ffsim.sample_state_vector(psi, norb=norb, nelec=nelec, shots=20000, seed=rng,
                                          bitstring_type=ffsim.BitstringType.INT)
    assert np.array_equal(samples, samples_b)
    ci, _ = lin_ci_strings(samples_b, norb, nelec, copy.deepcopy(rng), shots=5000, spb=spb, n_batches=nb,
                           max_dim=md_, subset_rng=np.random.default_rng(12345))
    es = [solve_sci(c, ham.one_body_tensor, ham.two_body_tensor, norb=norb, nelec=nelec, spin_sq=0.0).energy
          for c in ci]
    e_a = [e for e, _ in hist[0]]
    print("her one-call batch energies:", np.round(e_a, 10), "shapes", [s for _, s in hist[0]])
    print("staged batch energies      :", np.round(es, 10), "dims", [len(c[0]) for c in ci])
    print(f"max |dE| = {np.max(np.abs(np.array(e_a) - np.array(es))):.2e};  returned min {res.energy:.10f}"
          f" vs staged min {min(es):.10f}")


# ===================================================================== summary
def cmd_report(a):
    st = {}
    for ln in open(RES / "lucj_states.jsonl"):
        r = json.loads(ln)
        st[(r["name"], r["cand"], r["seed"])] = r
    diag = {}
    for f in (WORK / "diag").glob("*.json"):
        name, cand, s, p = f.stem.split("__")
        diag.setdefault((name, cand, int(s[1:])), {})[p] = json.load(open(f))["batches"]
    rows = []
    for k in sorted(st):
        r = dict(st[k])
        d = diag.get(k, {})
        for p in PROTOCOLS:
            if p in d:
                e = np.array([x["energy"] for x in d[p]])
                r[p] = dict(e=e.tolist(), mean=float(e.mean()), min=float(e.min()), max=float(e.max()),
                            dims=[x["dim_a"] for x in d[p]], s2=[x["spin_square"] for x in d[p]],
                            t_solve=[x["t_solve"] for x in d[p] if not x.get("dedup")],
                            n_distinct=len({x["sha1"] for x in d[p]}))
        rows.append(r)
    fci = json.load(open(RES / "fci_refs.json")) if (RES / "fci_refs.json").exists() else {}
    out = dict(sqd=SQD, protocols=PROTOCOLS, rows=rows, fci=fci)
    json.dump(out, open(RES / "summary.json", "w"), indent=1)
    print(f"-> {RES / 'summary.json'}  ({len(rows)} rows)")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("candidates")
    c.add_argument("--names-file", required=True)
    c.add_argument("--labels-dir", required=True)
    c.add_argument("--tag", required=True)
    c.add_argument("--with-cdf", action="store_true")
    c.add_argument("--threads", type=int, default=8)
    s = sub.add_parser("sample")
    s.add_argument("--cands", required=True)
    s.add_argument("--only-names", nargs="*")
    s.add_argument("--cands-only", nargs="*")
    s.add_argument("--seeds", type=int, nargs="+", default=[0])
    s.add_argument("--max-mem-gb", type=float, default=9.0)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--ci-workers", type=int, default=5)
    d = sub.add_parser("diag")
    d.add_argument("--threads", type=int, default=8)
    d.add_argument("--only-names", nargs="*")
    d.add_argument("--only-cands", nargs="*")
    d.add_argument("--skip-cands", nargs="*")
    d.add_argument("--names-file", default=None, help="restrict to these molecules (adds to --only-names)")
    d.add_argument("--seeds", type=int, nargs="*")
    d.add_argument("--protocols", nargs="*")
    d.add_argument("--stop-file", default=None)
    d.add_argument("--idle-exit-s", type=float, default=0, help="poll for new sample files until idle this long")
    f = sub.add_parser("fci")
    f.add_argument("--names-file", default="pretrain/rl/small_val.txt")
    f.add_argument("--names", nargs="*")
    f.add_argument("--threads", type=int, default=16)
    f.add_argument("--max-memory-mb", type=int, default=60000)
    f.add_argument("--conv-tol", type=float, default=1e-10)
    f.add_argument("--redo", action="store_true")
    f.add_argument("--verbose", action="store_true")
    t = sub.add_parser("selftest")
    t.add_argument("--threads", type=int, default=4)
    t.add_argument("--gpu-check", action="store_true")
    sub.add_parser("report")
    a = ap.parse_args()
    {"candidates": cmd_candidates, "sample": cmd_sample, "diag": cmd_diag, "fci": cmd_fci,
     "selftest": cmd_selftest, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
