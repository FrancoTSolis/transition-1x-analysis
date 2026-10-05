#!/usr/bin/env python3
"""Task 2: per-molecule derivative-free optimization of the LUCJ parameters, 500 objective evaluations.

  objective  lucj : exact LUCJ energy (pretrain.rl.gpu_energy, complex64)                      -> RL upper bound
             qsci : QSCI energy of 10^4 exact samples (Lin et al.'s TN-optimization objective, her settings,
                    samples drawn from the exact state instead of her chi=50 MPS, fixed sampling seed as her
                    quimb sample(seed=0); subsampling rng persistent across evaluations as hers)
  optimizer  nomad: PyNomad 4.5.1 PSD-MADS exactly as her lucj_compressed_t2_nomad_r24.py / the paper
                    (MAX_BB_EVAL 500, 20 variables and 20 evaluations per subproblem, 4 subproblems, no bounds,
                    default mesh); evaluations are serialized (one GPU)
             spsa : first-order SPSA (Spall gains alpha 0.602, gamma 0.101, A = 10% of iterations), learning
                    rate calibrated on n_cal gradient samples at x0 (they count toward the budget and their mean
                    is the first step), the best evaluated point is returned
  start      label (optimize=True, her start) | rl4n29f (network + RL, one call)

  CUDA_VISIBLE_DEVICES=3 python -m pretrain.followups.task2_opt --name C2H3N_rxn2858_P --start label \
      --objective lucj --optimizer nomad [--budget 500] [--worker-threads 2]

Outputs: <results>/<objective>/<name>__<start>__<optimizer>[__tag].{json,npz}; per-evaluation history
<logs>/<objective>/<same>.jsonl (written as it runs; best x checkpointed to the npz on every improvement).
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"          # before anything loads an OpenMP runtime (PyNomad PSD-MADS)
for _v in ("MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import threading  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pretrain.followups import task2_common as C  # noqa: E402
from pretrain.followups.task2_worker import _read, _write  # noqa: E402


class Remote:
    def __init__(self, name, threads, max_mem_gb, entropy, log_path, sci=True):
        env = dict(os.environ)
        for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS",
                  "NUMBA_NUM_THREADS"):
            env[v] = str(threads)
        self.log = open(log_path, "a")
        self.p = subprocess.Popen(
            [sys.executable, "-m", "pretrain.followups.task2_worker", "--name", name, "--threads", str(threads),
             "--max-mem-gb", str(max_mem_gb), "--entropy", str(entropy)] + ([] if sci else ["--no-sci"]),
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log, cwd=str(C.ROOT), env=env)
        self.info = _read(self.p.stdout)
        if not self.info or not self.info.get("ready"):
            raise RuntimeError(f"worker failed to start: {self.info}")

    def call(self, req):
        _write(self.p.stdin, req)
        r = _read(self.p.stdout)
        if r is None:
            raise RuntimeError("worker died")
        if "error" in r:
            raise RuntimeError(r["error"] + "\n" + r.get("tb", ""))
        return r

    def close(self):
        try:
            _write(self.p.stdin, {"cmd": "quit"})
            self.p.wait(timeout=60)
        except Exception:  # noqa: BLE001
            self.p.kill()


class StopOpt(Exception):
    """Wall-time limit reached or stop file present: no further evaluations."""


class Tracker:
    """Counts evaluations, keeps the history and the best point (checkpointed).  Stops (StopOpt) once the wall
    time exceeds --max-wall-h or the file <logs>/<objective>/<stem>.stop exists; the evaluations done so far
    are kept and reported (stopped_early)."""

    def __init__(self, remote, args, hist_path: Path, npz_path: Path, x0):
        self.r, self.a = remote, args
        self.n, self.best_f, self.best_x, self.best_n = 0, np.inf, None, -1
        self.fs, self.t0 = [], time.time()
        self.hist_path, self.npz_path, self.x0 = hist_path, npz_path, np.asarray(x0)
        hist_path.parent.mkdir(parents=True, exist_ok=True)
        self.hf = open(hist_path, "w")
        self.settings = dict(C.QSCI_OPT, shots=args.shots)
        self.stop_file = hist_path.with_suffix(".stop")
        self.stopped = None

    def __call__(self, x):
        if self.stopped is None:
            if self.a.max_wall_h and time.time() - self.t0 > 3600 * self.a.max_wall_h:
                self.stopped = f"wall time > {self.a.max_wall_h} h after {self.n} evaluations"
            elif self.stop_file.exists():
                self.stopped = f"stop file after {self.n} evaluations"
        if self.stopped is not None:
            raise StopOpt(self.stopped)
        x = np.asarray(x, dtype=np.float64)
        if self.a.objective == "lucj":
            r = self.r.call({"cmd": "lucj", "x": x})
        else:
            r = self.r.call({"cmd": "qsci", "x": x, "shots": self.a.shots, "sample_seed": self.a.sample_seed,
                             "settings": self.settings, "verify": True})
            if abs(r.get("d_full", 0.0)) > 1e-5:
                print(f"WARNING eval {self.n + 1}: SCI energy {r['f']:.8f} vs full-space Rayleigh quotient "
                      f"{r['E_full']:.8f} (asym {r.get('asym')})", flush=True)
        f = float(r["f"])
        self.n += 1
        self.fs.append(f)
        rec = {"n": self.n, "f": f, "t_eval": r["t"], "wall": time.time() - self.t0,
               "dx": float(np.linalg.norm(x - self.x0))}
        if self.a.objective == "qsci":
            rec.update(E_mean=r["E_mean"], dim=r["dims"][0], n_unique=r["n_unique"], d_full=r.get("d_full"),
                       asym=r.get("asym"), t_sci=r.get("t_sci"))
        if f < self.best_f:
            self.best_f, self.best_x, self.best_n = f, x.copy(), self.n
            tmp = self.npz_path.with_suffix(".tmp.npz")
            np.savez(tmp, x_best=self.best_x, x0=self.x0, f_best=f, n=self.n)
            tmp.replace(self.npz_path)
        rec["best"] = self.best_f
        self.hf.write(json.dumps(rec) + "\n")
        self.hf.flush()
        return f


def run_nomad(tr: Tracker, x0, params):
    import PyNomad
    lock = threading.Lock()

    def bb(p):
        try:
            x = np.array([p.get_coord(i) for i in range(p.size())])
            with lock:
                f = tr(x)
            p.setBBO(str(f).encode("UTF-8"))
            return 1
        except StopOpt:
            return 0                                  # failed evaluation: NOMAD runs out its budget at once
        except Exception:  # noqa: BLE001
            traceback.print_exc()
            return 0

    res = PyNomad.optimize(bb, list(map(float, x0)), [], [], params)
    return {k: (v if k != "x_best" else None) for k, v in res.items()}


def run_spsa(tr: Tracker, x0, budget, c, target, n_cal, seed, alpha=0.602, gamma=0.101, a_fixed=None):
    info = {}
    try:
        return _spsa(tr, x0, budget, c, target, n_cal, seed, alpha, gamma, info, a_fixed)
    except StopOpt as e:
        info["stopped"] = str(e)
        return info


def _spsa(tr, x0, budget, c, target, n_cal, seed, alpha, gamma, info, a_fixed=None):
    rng = np.random.default_rng(seed)
    n = len(x0)
    x = np.asarray(x0, dtype=np.float64).copy()
    tr(x)                                             # f(x0)
    # calibration (counts toward the budget): gradient samples at x0
    gs, mags = [], []
    for _ in range(n_cal):
        d = rng.choice([-1.0, 1.0], size=n)
        fp, fm = tr(x + c * d), tr(x - c * d)
        gs.append((fp - fm) / (2 * c) * d)
        mags.append(abs(fp - fm) / (2 * c))
    n_iter = (budget - 1 - 2 * n_cal - 1) // 2        # leave one evaluation for the final iterate
    A = 0.1 * n_iter
    a = target * (A + 1) ** alpha / max(np.mean(mags), 1e-12)
    if a_fixed:                                       # gain given (e.g. calibrated at the label start)
        info["a_calibrated"] = float(a)
        a = float(a_fixed)
    info.update(a=float(a), A=A, c=c, target=target, n_cal=n_cal, n_iter=n_iter, mean_grad_mag=float(np.mean(mags)))
    x = x - a / (A + 1) ** alpha * np.mean(gs, axis=0)  # first step from the calibration gradients
    for k in range(1, n_iter + 1):
        ak, ck = a / (k + 1 + A) ** alpha, c / (k + 1) ** gamma
        d = rng.choice([-1.0, 1.0], size=n)
        fp, fm = tr(x + ck * d), tr(x - ck * d)
        x = x - ak * (fp - fm) / (2 * ck) * d
    tr(x)                                             # final iterate
    info["x_final_f"] = tr.fs[-1]
    return info


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--start", default="label")
    ap.add_argument("--objective", choices=["lucj", "qsci"], default="lucj")
    ap.add_argument("--optimizer", choices=["nomad", "spsa"], default="nomad")
    ap.add_argument("--budget", type=int, default=500)
    ap.add_argument("--shots", type=int, default=C.QSCI_OPT["shots"])
    ap.add_argument("--sample-seed", type=int, default=0)
    ap.add_argument("--entropy", type=int, default=0)
    ap.add_argument("--spsa-c", type=float, default=0.005)
    ap.add_argument("--spsa-target", type=float, default=0.002, help="per-coordinate size of the first step")
    ap.add_argument("--spsa-ncal", type=int, default=10)
    ap.add_argument("--spsa-a", type=float, default=0.0, help="fixed SPSA gain a (0: calibrate from the target)")
    ap.add_argument("--spsa-a-from", default=None,
                    help="stem of a finished SPSA result (same objective) whose calibrated gain a is used, e.g. "
                         "<name>__label__spsa: the step-size scale of the label start applied to another start")
    ap.add_argument("--seed", type=int, default=0, help="SPSA perturbation seed")
    ap.add_argument("--nomad-extra", action="append", default=[], help="extra NOMAD parameter lines")
    ap.add_argument("--worker-threads", type=int, default=2)
    ap.add_argument("--max-mem-gb", type=float, default=5.0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--max-wall-h", type=float, default=0.0, help="stop evaluating after this many hours (0: none)")
    ap.add_argument("--finalize", action="store_true",
                    help="write the result JSON of a killed run from its checkpoint (npz + history)")
    args = ap.parse_args()
    if args.finalize:
        return finalize(args)

    stem = f"{args.name}__{args.start}__{args.optimizer}" + (f"__{args.tag}" if args.tag else "")
    res_dir, log_dir = C.RESULTS / args.objective, C.LOGS / args.objective
    res_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    U0, Z0 = C.start_uz(args.name, args.start)
    x0 = C.uz_to_x(U0, Z0)
    t_start = time.time()
    remote = Remote(args.name, args.worker_threads, args.max_mem_gb, args.entropy, log_dir / f"{stem}.worker.log",
                    sci=args.objective == "qsci")
    info = remote.info
    tr = Tracker(remote, args, log_dir / f"{stem}.jsonl", res_dir / f"{stem}.npz", x0)
    out = {"name": args.name, "start": args.start, "objective": args.objective, "optimizer": args.optimizer,
           "budget": args.budget, "args": vars(args), "norb": info["norb"], "nelec": info["nelec"],
           "e_hf": info["e_hf"], "e_ccsd": info["e_ccsd"], "n_params": len(x0), "gpu": info["gpu"],
           "host": os.uname().nodename, "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    if args.objective == "qsci":
        out["qsci_settings"] = dict(C.QSCI_OPT, shots=args.shots)
    try:
        if args.optimizer == "nomad":
            params = list(C.NOMAD_PARAMS)
            params = [p if not p.startswith("MAX_BB_EVAL") else f"MAX_BB_EVAL {args.budget}" for p in params]
            params += args.nomad_extra
            out["nomad_params"] = params
            out["nomad"] = run_nomad(tr, x0, params)
        else:
            a_fixed = args.spsa_a or None
            if args.spsa_a_from:
                src = json.load(open(res_dir / f"{args.spsa_a_from}.json"))
                if src["budget"] != args.budget:
                    raise ValueError("--spsa-a-from: the source run must have the same budget (same A)")
                a_fixed = float(src["spsa"]["a"])
                out["spsa_a_from"] = args.spsa_a_from
            out["spsa"] = run_spsa(tr, x0, args.budget, args.spsa_c, args.spsa_target, args.spsa_ncal, args.seed,
                                   a_fixed=a_fixed)
        t_opt = time.time() - t_start
        if tr.stopped:
            out["stopped_early"] = tr.stopped
        # variational (exact LUCJ) energy of the start and of the best point; not counted in the budget
        E0 = remote.call({"cmd": "lucj", "x": x0})["f"]
        Eb = remote.call({"cmd": "lucj", "x": tr.best_x})["f"]
    finally:
        remote.close()
    eh, ec = info["e_hf"], info["e_ccsd"]
    fs = np.array(tr.fs)
    out.update(n_evals=tr.n, wall_s=t_opt, t_eval_mean=None,
               f0=float(fs[0]), f_best=float(tr.best_f), best_at=tr.best_n,
               E_var0=float(E0), E_var_best=float(Eb), corr_var0=C.corr_pct(E0, eh, ec),
               corr_var_best=C.corr_pct(Eb, eh, ec), corr_f0=C.corr_pct(fs[0], eh, ec),
               corr_f_best=C.corr_pct(tr.best_f, eh, ec),
               best_trace=[float(v) for v in np.minimum.accumulate(fs)],
               dx_best=float(np.linalg.norm(tr.best_x - x0)), finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    hist = [json.loads(ln) for ln in open(log_dir / f"{stem}.jsonl")]
    out["t_eval_mean"] = float(np.mean([h["t_eval"] for h in hist]))
    C.dumpj(out, res_dir / f"{stem}.json")
    print(f"{stem}: {tr.n} evals, {t_opt/60:.1f} min; f {fs[0]:.6f} -> {tr.best_f:.6f} "
          f"({out['corr_f0']:.2f} -> {out['corr_f_best']:.2f} %corr); var {out['corr_var0']:.2f} -> "
          f"{out['corr_var_best']:.2f}", flush=True)


def finalize(args):
    """Result JSON of a run that was killed (or is to be cut short): from the history (.jsonl) and the best-point
    checkpoint (.npz); exact variational energies of x0 and x_best from a fresh worker.  Marked truncated."""
    stem = f"{args.name}__{args.start}__{args.optimizer}" + (f"__{args.tag}" if args.tag else "")
    res_dir, log_dir = C.RESULTS / args.objective, C.LOGS / args.objective
    hist = [json.loads(ln) for ln in open(log_dir / f"{stem}.jsonl")]
    z = np.load(res_dir / f"{stem}.npz")
    x0, xb = z["x0"], z["x_best"]
    remote = Remote(args.name, args.worker_threads, args.max_mem_gb, args.entropy, log_dir / f"{stem}.worker.log",
                    sci=False)
    try:
        E0 = remote.call({"cmd": "lucj", "x": x0})["f"]
        Eb = remote.call({"cmd": "lucj", "x": xb})["f"]
    finally:
        remote.close()
    info = remote.info
    eh, ec = info["e_hf"], info["e_ccsd"]
    fs = np.array([h["f"] for h in hist])
    out = {"name": args.name, "start": args.start, "objective": args.objective, "optimizer": args.optimizer,
           "budget": args.budget, "args": vars(args), "norb": info["norb"], "nelec": info["nelec"],
           "e_hf": eh, "e_ccsd": ec, "n_params": len(x0), "gpu": info["gpu"], "host": os.uname().nodename,
           "truncated": f"finalized from checkpoint after {len(hist)} evaluations",
           "n_evals": len(hist), "wall_s": hist[-1]["wall"], "t_eval_mean": float(np.mean([h["t_eval"] for h in hist])),
           "f0": float(fs[0]), "f_best": float(fs.min()), "best_at": int(np.argmin(fs)) + 1,
           "E_var0": float(E0), "E_var_best": float(Eb), "corr_var0": C.corr_pct(E0, eh, ec),
           "corr_var_best": C.corr_pct(Eb, eh, ec), "corr_f0": C.corr_pct(fs[0], eh, ec),
           "corr_f_best": C.corr_pct(fs.min(), eh, ec), "best_trace": [float(v) for v in np.minimum.accumulate(fs)],
           "dx_best": float(np.linalg.norm(xb - x0)), "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
    if abs(float(z["f_best"]) - fs.min()) > 1e-9:
        out["warning"] = f"checkpoint f_best {float(z['f_best'])} != history min {fs.min()}"
    if args.objective == "qsci":
        out["qsci_settings"] = dict(C.QSCI_OPT, shots=args.shots)
    C.dumpj(out, res_dir / f"{stem}.json")
    print(f"{stem}: finalized ({len(hist)} evals); f {fs[0]:.6f} -> {fs.min():.6f}; var {out['corr_var0']:.2f} -> "
          f"{out['corr_var_best']:.2f}", flush=True)


if __name__ == "__main__":
    main()
