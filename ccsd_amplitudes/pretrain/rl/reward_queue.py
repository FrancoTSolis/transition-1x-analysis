#!/usr/bin/env python3
"""Shared-filesystem reward queue: GRPO drivers submit LUCJ energy tasks, GPU workers on any machine claim them.

Layout under <root> (on the shared /xuanwu-tank filesystem, visible from scai1-7):
    pending/n<norb>_<na>_<nb>__<task_id>.pkl      task waiting for a worker (size in the name for routing)
    claimed/<worker>/<same file name>             task being computed (atomic rename = claim)
    done/<task_id>.pkl                            result {"id", "E", "info", "worker", "t"}
    workers/<worker>.json                         heartbeat {"host", "gpu", "max_norb", "time", "n_done"}
Task pickle: {"id", "name", "U", "Z", "t1", "norb", "nelec", "kind"}.

Driver side:   ids = submit(root, tasks);  res = collect(root, ids)   (requeues tasks of dead workers)
Worker side:   python3 -m pretrain.rl.reward_queue worker --root R --worker-id W --max-norb 18 [--dtype complex64]

GPU-limit guard (user rule): on scai3-scai7 the user account may hold at most 2 GPUs per machine, counting ALL of
the user's processes (other projects included).  `python3 -m pretrain.rl.reward_queue guard` exits non-zero if
starting one more GPU process here would exceed the limit; launch scripts call it first.
"""
from __future__ import annotations

import argparse
import getpass
import json
import os
import pickle
import socket
import subprocess
import sys
import time
import uuid
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GPU_LIMIT_HOSTS = {"scai3", "scai4", "scai5", "scai6", "scai7"}
GPU_LIMIT = 2


def _dirs(root: Path):
    for d in ("pending", "claimed", "done", "workers"):
        (root / d).mkdir(parents=True, exist_ok=True)


def submit(root, tasks):
    """tasks: list of dicts with name, U, Z, t1, norb, nelec (tuple), kind ('exact' | 'tn'), optional opts.
    Returns the task ids in order."""
    root = Path(root)
    _dirs(root)
    ids = []
    for t in tasks:
        tid = uuid.uuid4().hex[:16]
        t = dict(t, id=tid)
        na, nb = t["nelec"]
        fn = f"n{t['norb']:02d}_{na}_{nb}__{tid}.pkl"
        tmp = root / "pending" / (fn + ".tmp")
        with open(tmp, "wb") as f:
            pickle.dump(t, f)
        os.replace(tmp, root / "pending" / fn)
        ids.append(tid)
    return ids


def _requeue_stale(root: Path, stale_s: float):
    """Move tasks claimed by workers whose heartbeat is older than stale_s back to pending."""
    now = time.time()
    for wdir in (root / "claimed").iterdir():
        hb = root / "workers" / f"{wdir.name}.json"
        try:
            alive = hb.exists() and now - json.loads(hb.read_text())["time"] < stale_s
        except Exception:  # noqa: BLE001
            alive = False
        if alive:
            continue
        for f in wdir.glob("*.pkl"):
            try:
                os.replace(f, root / "pending" / f.name)
            except FileNotFoundError:
                pass


def collect(root, ids, timeout=None, poll=2.0, stale_s=300.0, verbose=False):
    """Wait for all ids; returns {id: (E, info)}.  E = nan on worker error."""
    root = Path(root)
    want, out = set(ids), {}
    t0, last_requeue = time.time(), 0.0
    while want:
        for tid in list(want):
            f = root / "done" / f"{tid}.pkl"
            if f.exists():
                try:
                    with open(f, "rb") as fh:
                        r = pickle.load(fh)
                except (EOFError, pickle.UnpicklingError):
                    continue                                   # still being written
                out[tid] = (r["E"], r.get("info", {}))
                want.discard(tid)
                f.unlink(missing_ok=True)
        if not want:
            break
        if time.time() - last_requeue > 30:
            _requeue_stale(root, stale_s)
            last_requeue = time.time()
        if timeout is not None and time.time() - t0 > timeout:
            for tid in want:
                out[tid] = (float("nan"), {"error": "timeout"})
            break
        if verbose and int(time.time() - t0) % 60 == 0:
            print(f"  [queue] waiting for {len(want)} results", flush=True)
        time.sleep(poll)
    return out


def user_gpus(user=None):
    """Set of GPU indices on this machine with a compute process owned by `user` (None on query failure)."""
    user = user or getpass.getuser()
    try:
        m = subprocess.run(["nvidia-smi", "--query-gpu=index,uuid", "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=30).stdout
        uuid2idx = {ln.split(",")[1].strip(): ln.split(",")[0].strip() for ln in m.strip().splitlines()}
        a = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,gpu_uuid", "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=30).stdout
    except Exception:  # noqa: BLE001
        return None
    gpus = set()
    for ln in a.strip().splitlines():
        if not ln.strip():
            continue
        pid, g = [x.strip() for x in ln.split(",")]
        try:
            owner = subprocess.run(["ps", "-o", "user=", "-p", pid], capture_output=True, text=True).stdout.strip()
        except Exception:  # noqa: BLE001
            owner = ""
        if owner == user:
            gpus.add(uuid2idx.get(g, g))
    return gpus


def user_gpu_count(user=None):
    g = user_gpus(user)
    return 99 if g is None else len(g)


def guard(extra=1, gpu=None):
    """Would one more process (on GPU index `gpu`, if given) keep the user within GPU_LIMIT GPUs on this host?"""
    host = socket.gethostname().split(".")[0]
    if host not in GPU_LIMIT_HOSTS:
        return True, host, None
    g = user_gpus()
    if g is None:
        return False, host, None
    if gpu is not None and str(gpu) in g:
        extra = 0                                       # that GPU is already counted
    return len(g) + extra <= GPU_LIMIT, host, len(g)


class EngineCache:
    """LRU cache of per-molecule energy engines (Hamiltonian tensors + string tables live on the GPU)."""

    def __init__(self, factory, max_items=4):
        self.factory, self.max_items, self.d = factory, max_items, OrderedDict()

    def get(self, name, norb, nelec):
        if name in self.d:
            self.d.move_to_end(name)
            return self.d[name]
        while len(self.d) >= self.max_items:
            _, old = self.d.popitem(last=False)
            del old
            import gc
            gc.collect()
            try:
                import torch
                torch.cuda.empty_cache()
            except Exception:  # noqa: BLE001
                pass
        eng = self.factory(name, norb, nelec)
        self.d[name] = eng
        return eng


def default_exact_factory(dtype_name="complex64", max_mem_gb=None):
    def make(name, norb, nelec):
        import numpy as np
        import torch
        from pretrain.rl.gpu_energy import LUCJEnergyGPU
        d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
        return LUCJEnergyGPU(d["one_body"], d["two_body"], float(d["constant"]), int(d["norb"]),
                             (int(d["nelec_a"]), int(d["nelec_b"])), device="cuda",
                             dtype=getattr(torch, dtype_name), max_mem_gb=max_mem_gb)
    return make


def tn_factory(chi=256, basis_cache=None, stack_mem_gb=2.0, zip_margin=1.5, impl="current", block2_threads=1):
    """Tensor-network (MPS) engines on this GPU: pretrain.rl.tn_energy.LUCJEnergyTN, complex64, block2 <H> on CPU.
    zip_margin (zip-up bond = margin * chi) is an accuracy knob like chi: 3.0 at chi 256 costs ~3x margin 1.5."""
    def make(name, norb, nelec):
        import tempfile
        import numpy as np
        if impl == "v1":          # frozen zip-up engine of the Oct-3 evaluations (pretrain/rl/tn_energy_v1.py)
            from pretrain.rl.tn_energy_v1 import LUCJEnergyTN
        else:
            from pretrain.rl.tn_energy import LUCJEnergyTN
        d = np.load(ROOT / "rhf_hamiltonians" / f"{name}.npz")
        return LUCJEnergyTN(d["one_body"], d["two_body"], float(d["constant"]), int(d["norb"]),
                            (int(d["nelec_a"]), int(d["nelec_b"])), max_bond=chi, device="cuda", name=name,
                            block2_threads=block2_threads, scratch=tempfile.mkdtemp(prefix="tnq_"),
                            stack_mem=int(stack_mem_gb * (1 << 30)), basis_cache=basis_cache,
                            zip_margin=zip_margin)
    return make


def _write_hb(path: Path, d: dict):
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(d))
    os.replace(tmp, path)                        # readers never see a half-written file


def _hb_process(path, info, n_done, stop, parent):
    while not stop.wait(20):
        if os.getppid() != parent:               # worker gone (even after SIGKILL): stop beating
            return
        try:
            _write_hb(path, {**info, "time": time.time(), "n_done": n_done.value})
        except OSError:
            pass


def worker(root, worker_id, max_norb, factory, kinds=("exact",), poll=1.0, max_items=4, idle_exit=None):
    import multiprocessing as mp
    import numpy as np
    root = Path(root)
    _dirs(root)
    mine = root / "claimed" / worker_id
    mine.mkdir(parents=True, exist_ok=True)
    cache = EngineCache(factory, max_items=max_items)
    host = socket.gethostname().split(".")[0]
    n_done, last_work = 0, time.time()
    # Heartbeat from a child process: one task can run far longer than collect()'s stale_s (TN energies at chi 256
    # take 5-30 min), and a busy worker must not look dead, or its task is requeued and computed twice.  A thread
    # is not enough: C extensions (block2's expectation value) hold the GIL for minutes.  Forked before any
    # engine (torch / CUDA) exists.
    hb_path = root / "workers" / f"{worker_id}.json"
    info = {"host": host, "gpu": os.environ.get("CUDA_VISIBLE_DEVICES"), "max_norb": max_norb, "pid": os.getpid()}
    ctx = mp.get_context("fork")
    n_done_v, stop_hb = ctx.Value("i", 0), ctx.Event()
    _write_hb(hb_path, {**info, "time": time.time(), "n_done": 0})
    hb_proc = ctx.Process(target=_hb_process, args=(hb_path, info, n_done_v, stop_hb, os.getpid()), daemon=True)
    hb_proc.start()

    def stop_heartbeat():
        stop_hb.set()
        hb_proc.join(timeout=30)

    last_guard = time.time()
    while True:
        # keep honouring the per-machine GPU limit while running: the user's other jobs may start later
        if host in GPU_LIMIT_HOSTS and time.time() - last_guard > 60:
            last_guard = time.time()
            n = user_gpu_count()
            if n > GPU_LIMIT:
                print(f"[{worker_id}] user now holds {n} GPUs on {host} (limit {GPU_LIMIT}): exiting", flush=True)
                stop_heartbeat()
                hb_path.unlink(missing_ok=True)
                break
        cands = sorted((root / "pending").glob("n*.pkl"))
        # largest systems first (they are the bottleneck), only sizes this GPU can hold
        cands = [f for f in cands if int(f.name[1:3]) <= max_norb]
        cands.sort(key=lambda f: -int(f.name[1:3]))
        claimed = None
        for f in cands:
            try:
                os.replace(f, mine / f.name)
                claimed = mine / f.name
                break
            except FileNotFoundError:
                continue                                       # another worker was faster
        if claimed is None:
            if idle_exit and time.time() - last_work > idle_exit:
                stop_heartbeat()
                hb_path.unlink(missing_ok=True)
                break
            time.sleep(poll)
            continue
        last_work = time.time()
        with open(claimed, "rb") as fh:
            t = pickle.load(fh)
        t0 = time.time()
        try:
            if t.get("kind", "exact") not in kinds:
                raise RuntimeError(f"worker kinds {kinds} cannot run {t.get('kind')}")
            eng = cache.get(t["name"], t["norb"], tuple(t["nelec"]))
            U = np.asarray(t["U"], dtype=np.complex128)
            W, _, Vh = np.linalg.svd(U)
            out = eng.energy(W @ Vh, np.asarray(t["Z"], dtype=np.float64), t1=t.get("t1"))
            extra = {}
            if isinstance(out, tuple):                 # TN engines return (E, info)
                out, tinfo = out
                extra = {k: float(v) for k, v in tinfo.items() if k in ("discarded_sum", "max_bond", "t_total")}
            E = float(out)
            info = {"t": time.time() - t0, "worker": worker_id, "host": host, **extra}
        except Exception as e:  # noqa: BLE001
            E, info = float("nan"), {"error": f"{type(e).__name__}: {e}", "worker": worker_id}
        eng = None        # only the cache may keep an engine alive (else eviction cannot free its GPU memory)
        tmp = root / "done" / f"{t['id']}.pkl.tmp"
        with open(tmp, "wb") as fh:
            pickle.dump({"id": t["id"], "E": E, "info": info, "worker": worker_id, "t": time.time() - t0}, fh)
        os.replace(tmp, root / "done" / f"{t['id']}.pkl")
        claimed.unlink(missing_ok=True)
        n_done += 1
        n_done_v.value = n_done


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("guard")
    g.add_argument("--extra", type=int, default=1)
    g.add_argument("--gpu", default=None, help="target GPU index (not counted twice if already in use by the user)")
    w = sub.add_parser("worker")
    w.add_argument("--root", required=True)
    w.add_argument("--worker-id", required=True)
    w.add_argument("--max-norb", type=int, default=16)
    w.add_argument("--dtype", default="complex64")
    w.add_argument("--max-mem-gb", type=float, default=None)
    w.add_argument("--max-items", type=int, default=1)
    w.add_argument("--idle-exit", type=float, default=None, help="exit after this many idle seconds")
    w.add_argument("--no-guard", action="store_true")
    w.add_argument("--kind", default="exact", choices=["exact", "tn"])
    w.add_argument("--chi", type=int, default=256, help="MPS bond dimension for --kind tn")
    w.add_argument("--zip-margin", type=float, default=1.5, help="zip-up margin for --kind tn")
    w.add_argument("--block2-threads", type=int, default=1, help="CPU threads for block2 <H> (--kind tn)")
    w.add_argument("--tn-impl", default="current", choices=["current", "v1"],
                   help="v1 = frozen zip-up engine (tn_energy_v1.py) used for the Oct-3 n29 evaluations")
    w.add_argument("--basis-cache", default=str(ROOT / "rl_runs" / "tn_basis_cache"))
    args = ap.parse_args()
    if args.cmd == "guard":
        ok, host, n = guard(args.extra, args.gpu)
        print(f"host {host}: user GPUs in use {n}, limit {GPU_LIMIT if host in GPU_LIMIT_HOSTS else 'none'} -> "
              f"{'OK' if ok else 'REFUSE'}")
        sys.exit(0 if ok else 1)
    if args.cmd == "worker":
        if not args.no_guard:
            ok, host, n = guard(1, os.environ.get("CUDA_VISIBLE_DEVICES"))
            if not ok:
                print(f"refusing to start on {host}: user already holds {n} GPUs (limit {GPU_LIMIT})", flush=True)
                sys.exit(2)
        sys.path.insert(0, str(ROOT))
        fac = (tn_factory(args.chi, args.basis_cache, zip_margin=args.zip_margin, impl=args.tn_impl,
                          block2_threads=args.block2_threads)
               if args.kind == "tn"
               else default_exact_factory(args.dtype, args.max_mem_gb))
        worker(args.root, args.worker_id, args.max_norb, fac, kinds=(args.kind,), max_items=args.max_items,
               idle_exit=args.idle_exit)


if __name__ == "__main__":
    main()
