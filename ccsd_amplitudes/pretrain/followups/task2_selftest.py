#!/usr/bin/env python3
"""Self-checks of task2_common (run on one GPU):

1. x <-> (U, Z) round trip and GPU energies of the start points vs the stored exact energies
   (energy_small_baselines.json for the label, energy_smallval_rl4n29.json for rl4n29f).
2. The GPU CI-matrix sampler vs ffsim.sample_state_vector on a random norb-8 LUCJ state
   (state agreement, bitstring convention, total-variation distance of 2e5-sample histograms).
3. Cost of one QSCI evaluation at Lin et al.'s optimization settings (1e4 samples) and of one scoring
   evaluation (1e5 samples) on the label start.

  CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=6 python -m pretrain.followups.task2_selftest --names A B --qsci
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.followups import task2_common as C  # noqa: E402


def check_sampler():
    import ffsim
    import torch
    from pretrain.rl.energy import make_ucj_op
    from pretrain.rl.gpu_energy import LUCJEnergyGPU
    rng = np.random.default_rng(1)
    n, k = 8, 4
    h1 = rng.normal(size=(n, n)); h1 = h1 + h1.T
    h2 = rng.normal(size=(n, n, n, n)) * 0.1
    h2 = h2 + h2.transpose(1, 0, 2, 3); h2 = h2 + h2.transpose(0, 1, 3, 2); h2 = h2 + h2.transpose(2, 3, 0, 1)
    eng = LUCJEnergyGPU(h1, h2, 0.0, n, (k, k), max_mem_gb=1.0)
    U = np.stack([ffsim.random.random_unitary(n, seed=s) for s in (2, 3)])
    Z = rng.normal(size=(2, n, n)) * 0.5; Z = Z + Z.transpose(0, 2, 1)
    t1 = rng.normal(size=(k, n - k)) * 0.1
    x = C.uz_to_x(U, Z)
    U2, Z2 = C.x_to_uz(x, n)
    op = make_ucj_op(Z2, U2, "square", t1=t1)
    vec = ffsim.apply_unitary(ffsim.hartree_fock_state(n, (k, k)), op, norb=n, nelec=(k, k))
    psi = eng.state(U2, Z2, t1=t1).clone()
    v = psi.cpu().numpy().reshape(-1).astype(np.complex128)
    ph = np.vdot(v, vec) / abs(np.vdot(v, vec))
    out = dict(state_maxdiff=float(np.abs(v * ph - vec).max()))
    shots = 200_000
    mine = C.sample_ci_matrix(psi, n, k, shots, seed=5)
    ref = ffsim.sample_state_vector(vec, norb=n, nelec=(k, k), shots=shots, seed=7,
                                    bitstring_type=ffsim.BitstringType.INT)
    ref = np.asarray(ref, dtype=np.int64)
    keys = np.union1d(mine, ref)
    hm = np.searchsorted(keys, mine); hr = np.searchsorted(keys, ref)
    pm = np.bincount(hm, minlength=len(keys)) / shots
    pr = np.bincount(hr, minlength=len(keys)) / shots
    # exact probabilities in the same keys (ffsim index -> int via addresses_to_strings)
    strs = ffsim.addresses_to_strings(np.arange(len(vec)), norb=n, nelec=(k, k), bitstring_type=ffsim.BitstringType.INT)
    pex = dict(zip(np.asarray(strs, dtype=np.int64).tolist(), (np.abs(vec) ** 2).tolist()))
    pe = np.array([pex.get(int(s), 0.0) for s in keys])
    out.update(tv_mine_ref=float(0.5 * np.abs(pm - pr).sum()), tv_mine_exact=float(0.5 * np.abs(pm - pe).sum()),
               tv_ref_exact=float(0.5 * np.abs(pr - pe).sum()), n_keys=int(len(keys)),
               mine_outside_support=float(pm[pe == 0].sum()))
    eng.release()
    del eng
    torch.cuda.empty_cache()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True)
    ap.add_argument("--qsci", action="store_true")
    ap.add_argument("--score", action="store_true", help="also time one 1e5-sample scoring evaluation")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--max-mem-gb", type=float, default=5.0)
    ap.add_argument("--out", default=str(C.RESULTS / "selftest.json"))
    args = ap.parse_args()
    from pyscf import lib
    lib.num_threads(args.threads)
    import torch
    from pretrain.rl.gpu_energy import LUCJEnergyGPU
    res = {"sampler": check_sampler()}
    print(json.dumps(res["sampler"]), flush=True)
    stored = {}
    for f, tag in (("energy_small_baselines.json", "label"), ("energy_smallval_rl4n29.json", "rl4n29f")):
        pm = json.load(open(C.ROOT / "pretrain/opt_true/results" / f))["per_molecule"]
        for n, v in pm.items():
            if tag in v:
                stored[(n, tag)] = v[tag]["E"]
    for name in args.names:
        ham = C.load_ham(name)
        t1 = C.load_t1(name)
        eng = LUCJEnergyGPU.from_npz(name, max_mem_gb=args.max_mem_gb)
        r = {}
        for start in ("label", "rl4n29f"):
            U, Z = C.start_uz(name, start)
            t0 = time.time()
            E_direct = eng.energy(U, Z, t1=t1)
            t_e = time.time() - t0
            x = C.uz_to_x(U, Z)
            U2, Z2 = C.x_to_uz(x, ham["norb"])
            E_rt = eng.energy(U2, Z2, t1=t1)
            r[start] = dict(E_stored=stored.get((name, start)), E_gpu=E_direct, E_roundtrip=E_rt,
                            d_stored=None if (name, start) not in stored else E_direct - stored[(name, start)],
                            d_roundtrip=E_rt - E_direct, t_energy=t_e, n_params=len(x),
                            corr=C.corr_pct(E_direct, ham["e_hf"], ham["e_ccsd"]))
            print(name, start, json.dumps(r[start]), flush=True)
            if args.qsci and start == "label":
                for shots, key in ((10_000, "qsci_opt"), (100_000, "qsci_score")):
                    if shots == 100_000 and not args.score:
                        continue
                    t0 = time.time()
                    psi = eng.state(U2, Z2, t1=t1)
                    torch.cuda.synchronize()
                    t_state = time.time() - t0
                    t0 = time.time()
                    bs = C.sample_ci_matrix(psi, ham["norb"], ham["nelec"][0], shots, seed=0)
                    t_samp = time.time() - t0
                    q = C.qsci_energy(ham, bs, dict(C.QSCI_OPT, shots=shots), np.random.default_rng(0))
                    q.update(t_state=t_state, t_sample=t_samp, corr=C.corr_pct(q["E"], ham["e_hf"], ham["e_ccsd"]),
                             corr_mean=C.corr_pct(q["E_mean"], ham["e_hf"], ham["e_ccsd"]))
                    r[key] = q
                    print(name, key, json.dumps({k: v for k, v in q.items() if k != "E_batches"}), flush=True)
        res[name] = r
        eng.release()
        del eng
        torch.cuda.empty_cache()
    C.dumpj(res, Path(args.out))


if __name__ == "__main__":
    main()
