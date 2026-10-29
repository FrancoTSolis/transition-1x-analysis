"""Direct optimization of the exact LUCJ energy (GPU adjoint gradient) from a chosen initialization.

Inits (LUCJ, n_reps = 2, square connectivity, final orbital rotation = t1, fixed):
  zero        U = I, Z = 0  (the state is the t1-rotated HF determinant; a stationary point)
  random      log U ~ 0.1 N(0,1) (complex antihermitian), J ~ 0.1 N(0,1)   (--rand-seed)
  truncated   Lin et al.'s truncated t2 init = ffsim from_t_amplitudes, optimize=False
  cdf_lin     Lin et al.'s compressed DF (ffsim optimize=True, no regularization)
  label       our optimize=True label (regularized compressed DF)
  pre4        pretrained network (4 recycles)        rl4L, rl4n29f   network + energy RL
  (candidates from runs_ot/energy_tasks/task1_sqd_vs_lin/cands_smallval.pkl)

Methods:
  lbfgs   scipy L-BFGS-B (maxcor 20) with exact energy + gradient
  trn     Avijit's trust-region Newton driver (avijit_simulator/UCJ_optimizer.py, RUCJ.kernel): Steihaug-CG in
          a trust region, radius 0.25 (max 1.0), accept if rho >= 0.1, shrink x0.25 if rho < 0.25, grow x2
          if rho > 0.75 at the boundary; Hessian-vector products by forward differences of exact gradients.
Every energy/gradient evaluation is logged (JSONL); the budget is counted in gradient evaluations.

    CUDA_VISIBLE_DEVICES=2 python -m ucj_init.opt_exact --name C2H3N_rxn2858_TS --init rl4n29f --method lbfgs \
        --max-evals 300
"""
from __future__ import annotations

import argparse
import json
import pickle
import time
from pathlib import Path

import numpy as np
import scipy.optimize

from ucj_init.lucj_grad_gpu import LUCJGradGPU, n_params, x_from_uz

ROOT = Path(__file__).resolve().parents[1]
CANDS = ROOT / "runs_ot" / "energy_tasks" / "task1_sqd_vs_lin" / "cands_smallval.pkl"
RES = ROOT / "pretrain" / "opt_true" / "results" / "ucj_init" / "exact_opt"
LOGS = ROOT / "rl_runs" / "ucj_init" / "exact_opt"


def load_cands(name, path=CANDS):
    d = pickle.load(open(path, "rb"))
    out = {k[1]: (np.asarray(U), np.asarray(Z), np.asarray(t1)) for k, n, U, Z, t1 in d["tasks"] if n == name}
    if not out:
        raise KeyError(name)
    return out


def init_x(init, cands, norb, rand_seed=0, rand_scale=0.1):
    npar = n_params(norb, 2)
    if init == "zero":
        return np.zeros(npar)
    if init == "random":
        rng = np.random.default_rng(rand_seed)
        return rand_scale * rng.normal(size=npar)
    U, Z, _ = cands[init]
    return x_from_uz(U, Z)


class JastrowOnly:
    """Restrict an Objective to the Jastrow parameters (orbital rotations frozen at the start)."""

    def __init__(self, obj, x0, norb, n_reps=2):
        self.obj, self.x0 = obj, x0.copy()
        per = norb * norb + 2 * norb - 1
        self.idx = np.concatenate([np.arange(r * per + norb * norb, (r + 1) * per) for r in range(n_reps)])

    def full(self, y):
        x = self.x0.copy()
        x[self.idx] = y
        return x

    def __call__(self, y):
        E, g = self.obj(self.full(y))
        return E, g[self.idx]


class Objective:
    def __init__(self, eng, log_path, max_evals):
        self.eng, self.max_evals = eng, max_evals
        self.n = 0
        self.best = (np.inf, None)
        self.log = open(log_path, "w")
        self.t0 = time.time()
        self.cache = {}

    def __call__(self, x):
        key = x.tobytes()
        if key in self.cache:
            return self.cache[key]
        if self.n >= self.max_evals:
            raise StopIteration
        E, g = self.eng.energy_grad(x)
        self.n += 1
        if E < self.best[0]:
            self.best = (E, x.copy())
        rec = dict(i=self.n, E=E, pct=self.eng.corr_pct(E), gnorm=float(np.linalg.norm(g)),
                   ginf=float(np.abs(g).max()), t=time.time() - self.t0)
        self.log.write(json.dumps(rec) + "\n")
        self.log.flush()
        self.cache = {key: (E, g)}
        return E, g


def steihaug(gradient, hessp, radius, max_cycle):
    """avijit_simulator/UCJ_optimizer._steihaug (same tolerances)."""
    step = np.zeros_like(gradient)
    residual = gradient.copy()
    direction = -residual
    rr = residual @ residual
    tol = min(0.5, np.sqrt(np.linalg.norm(gradient))) * np.linalg.norm(gradient)
    if np.sqrt(rr) <= tol:
        return step, False, 0

    def boundary(s, d):
        aa, bb, cc = d @ d, 2 * s @ d, s @ s - radius ** 2
        return (-bb + np.sqrt(max(0.0, bb * bb - 4 * aa * cc))) / (2 * aa)
    for cycle in range(1, max_cycle + 1):
        hd = hessp(direction)
        curv = direction @ hd
        if curv <= 0:
            return step + boundary(step, direction) * direction, True, cycle
        alpha = rr / curv
        trial = step + alpha * direction
        if np.linalg.norm(trial) >= radius:
            return step + boundary(step, direction) * direction, True, cycle
        step = trial
        residual = residual + alpha * hd
        rr_new = residual @ residual
        if np.sqrt(rr_new) <= tol:
            return step, False, cycle
        direction = -residual + (rr_new / rr) * direction
        rr = rr_new
    return step, np.linalg.norm(step) >= 0.99 * radius, max_cycle


def trust_region_newton(obj, x0, inner=20, fd_h=1e-4, radius=0.25, radius_max=1.0, accept=0.1,
                        max_reject=10, gtol=1e-5):
    """RUCJ.kernel's macro-iteration logic on a generic objective (exact E, g; FD Hessian-vector)."""
    x = x0.copy()
    E, g = obj(x)
    hist = []
    while True:
        if np.abs(g).max() < gtol:
            return x, True, hist

        def hessp(v):
            nv = np.linalg.norm(v)
            if nv == 0:
                return np.zeros_like(v)
            _, gp = obj.eng.energy_grad(x + (fd_h / nv) * v)
            obj.n += 1
            if obj.n > obj.max_evals:
                raise StopIteration
            return (gp - g) * (nv / fd_h)
        accepted = False
        for attempt in range(max_reject):
            step, at_b, ncg = steihaug(g, hessp, radius, inner)
            hs = hessp(step)
            pred = -(g @ step + 0.5 * step @ hs)
            if pred <= 0 or not np.isfinite(pred):
                radius *= 0.25
                continue
            Et, gt = obj(x + step)
            actual = E - Et
            ratio = actual / pred
            if ratio < 0.25:
                radius *= 0.25
            elif ratio > 0.75 and at_b:
                radius = min(2 * radius, radius_max)
            if actual > 0 and ratio >= accept:
                x, E, g = x + step, Et, gt
                hist.append(dict(n=obj.n, E=E, radius=radius, ratio=ratio, ncg=ncg))
                accepted = True
                break
        if not accepted:
            return x, False, hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--init", required=True)
    ap.add_argument("--method", default="lbfgs", choices=["lbfgs", "trn"])
    ap.add_argument("--max-evals", type=int, default=300)
    ap.add_argument("--rand-seed", type=int, default=0)
    ap.add_argument("--rand-scale", type=float, default=0.1)
    ap.add_argument("--inner", type=int, default=20)
    ap.add_argument("--tag", default="")
    ap.add_argument("--cands", default=str(CANDS))
    ap.add_argument("--jastrow-only", action="store_true", help="freeze the orbital rotations at the start")
    args = ap.parse_args()
    cands = load_cands(args.name, args.cands)
    t1 = next(iter(cands.values()))[2]
    eng = LUCJGradGPU.from_npz(args.name, t1=t1)
    init = args.init + (f"{args.rand_seed}" if args.init == "random" else "")
    run = f"{args.name}__{init}__{args.method}{'J' if args.jastrow_only else ''}{args.tag}"
    RES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    x0 = init_x(args.init, cands, eng.norb, args.rand_seed, args.rand_scale)
    obj = Objective(eng, LOGS / f"{run}.jsonl", args.max_evals)
    t0 = time.time()
    status = "max_evals"
    E0, g0 = obj(x0)
    try:
        if args.jastrow_only:
            jo = JastrowOnly(obj, x0, eng.norb)
            res = scipy.optimize.minimize(jo, x0[jo.idx], jac=True, method="L-BFGS-B",
                                          options=dict(maxfun=args.max_evals, maxiter=10 * args.max_evals,
                                                       maxcor=20, gtol=1e-6, ftol=1e-13))
            status = f"scipy:{res.status}:{res.message}"
        elif args.method == "lbfgs":
            res = scipy.optimize.minimize(obj, x0, jac=True, method="L-BFGS-B",
                                          options=dict(maxfun=args.max_evals, maxiter=10 * args.max_evals,
                                                       maxcor=20, gtol=1e-6, ftol=1e-13))
            status = f"scipy:{res.status}:{res.message}"
        else:
            _, conv, _ = trust_region_newton(obj, x0, inner=args.inner)
            status = "converged" if conv else "tr_failed"
    except StopIteration:
        pass
    Eb, xb = obj.best
    out = dict(name=args.name, init=init, method=args.method + ("J" if args.jastrow_only else ""), n_evals=obj.n, seconds=time.time() - t0,
               status=status, e_hf=eng.e_hf, e_ccsd=eng.e_ccsd, E0=E0, pct0=eng.corr_pct(E0),
               gnorm0=float(np.linalg.norm(g0)), E_best=Eb, pct_best=eng.corr_pct(Eb), norb=eng.norb,
               n_params=eng.n_params, sec_per_eval=(time.time() - t0) / max(1, obj.n))
    np.savez(RES / f"{run}.npz", x0=x0, x_best=xb)
    json.dump(out, open(RES / f"{run}.json", "w"), indent=1)
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
