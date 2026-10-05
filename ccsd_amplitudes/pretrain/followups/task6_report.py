#!/usr/bin/env python3
"""Task 6 tables: the n29 GRPO run with dm chi-128 TN rewards (rl_runs/followups/task6_n29_rl_dm/run) next to the
Expanse v1 chi-64 run (rl_runs/grpo_n29_tn), and v1 chi-256 evaluations of its best / last policies on n29 val20 /
test24 next to the stored candidates of docs/rl_larger_molecules.md.

    python3 pretrain/followups/task6_report.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "pretrain/opt_true/results"
NEW = RES / "followups/task6_n29_rl_dm"
TD = ROOT / "rl_runs/followups/task6_n29_rl_dm"
RUN = TD / "run"
LAB = {"label": "optimize=True label", "pre4": "pretrained, 4 recycles", "rl4s": "small-set RL, 4 recycles",
       "rl4L": "start: norb 15-18 RL (step 25)", "rl4n29": "Expanse v1 chi-64 RL, step 15 (best)",
       "rl4n29f": "Expanse v1 chi-64 RL, step 30 (last)", "rl4n29dm": "this run: dm chi-128 RL, step 20 (best)",
       "rl4n29dmf": "this run: dm chi-128 RL, step 24 (last)"}


def load(p):
    p = Path(p)
    return json.load(open(p))["per_molecule"] if p.exists() else None


def merge(*pms):
    out = {}
    for pm in pms:
        for n, d in (pm or {}).items():
            out.setdefault(n, {}).update({c: v["corr_frac"] for c, v in d.items() if not v.get("error")})
    return out


def table(per, cands, ref="label"):
    names = sorted(n for n in per if ref in per[n])
    r = np.array([per[n][ref] for n in names])
    print(f"| candidate | n | mean % corr | median | min | vs {ref} (mean) | better than {ref} |")
    print("|:--|--:|--:|--:|--:|--:|--:|")
    for c in cands:
        ok = [i for i, n in enumerate(names) if c in per[n]]
        if not ok:
            continue
        v = np.array([per[names[i]][c] for i in ok])
        d = v - r[ok]
        cmp_ = "—" if c == ref else f"{100 * d.mean():+.1f}"
        bt = "—" if c == ref else f"{int((d > 0).sum())}/{len(ok)}"
        print(f"| {LAB.get(c, c)} | {len(ok)} | {100 * v.mean():.1f} | {100 * np.median(v):.1f} | "
              f"{100 * v.min():.1f} | {cmp_} | {bt} |")
    print()


def pairwise(per, a, b):
    names = sorted(n for n in per if a in per[n] and b in per[n])
    d = np.array([per[n][a] - per[n][b] for n in names])
    if len(d) == 0:
        return None
    return {"n": len(d), "mean_pp": 100 * d.mean(), "sem_pp": 100 * d.std(ddof=1) / np.sqrt(len(d)),
            "better": int((d > 0).sum())}


def curve(path, label):
    p = Path(path) / "log.jsonl"
    if not p.exists():
        return
    rows = [json.loads(ln) for ln in open(p)]
    # the Expanse log holds two failed starts: keep the last run (from the last step-0 val eval on)
    i0 = max(i for i, r in enumerate(rows) if r.get("type") == "eval_val" and r["step"] == 0)
    rows = rows[i0:]
    ev = {r["step"]: r for r in rows if r.get("type") == "eval_val"}
    tr = {r["step"]: r for r in rows if r.get("type") == "eval_train"}
    st = [r for r in rows if r.get("type") == "train"]
    steps = sorted(ev)
    print(f"**{label}**\n")
    print("| step | " + " | ".join(str(s) for s in steps) + " |")
    print("|:--|" + "--:|" * len(steps))
    print("| val (20) mean % corr | " + " | ".join(f"{100 * ev[s]['corr_mean']:.2f}" for s in steps) + " |")
    print("| val median | " + " | ".join(f"{100 * ev[s]['corr_median']:.2f}" for s in steps) + " |")
    print("| train (40) mean | " + " | ".join(f"{100 * tr[s]['corr_mean']:.2f}" if s in tr else "" for s in steps)
          + " |")
    if st:
        tr_ = np.array([r["t_reward"] for r in st])
        print(f"\n{len(st)} steps, reward time per step median {np.median(tr_):.0f} s (min {tr_.min():.0f}, max "
              f"{tr_.max():.0f}), failed energies {sum(r['n_fail'] for r in st)}; train reward mean per step: "
              + ", ".join(f"{100 * r['reward_mean']:.1f}" for r in st))
    t0, t1 = rows[0]["time"], rows[-1]["time"]
    print(f"wall clock (first eval to last record): {(t1 - t0) / 3600:.2f} h\n")


def drift(path, label):
    """Per-step train reward minus the start policy's (step-0 eval_train) value on the same molecules, and paired
    val / train changes against step 0.  The batches are replayed from the driver's numpy RNG (default_rng(seed),
    used only for rng.choice(train_names, batch_mols) once per step in grpo_slot.py)."""
    p = Path(path) / "log.jsonl"
    if not p.exists():
        return None
    rows = [json.loads(ln) for ln in open(p)]
    i0 = max(i for i, r in enumerate(rows) if r.get("type") == "eval_val" and r["step"] == 0)
    rows = rows[i0:]
    a = json.load(open(Path(path) / "args.json"))
    names = [ln.strip() for ln in open(ROOT / a["train_names"]) if ln.strip()]
    et = {r["step"]: r["per_mol"] for r in rows if r.get("type") == "eval_train"}
    ev = {r["step"]: r["per_mol"] for r in rows if r.get("type") == "eval_val"}
    tr = {r["step"]: r for r in rows if r.get("type") == "train"}
    rng = np.random.default_rng(a["seed"])
    dd = []
    for s in sorted(tr):
        b = list(rng.choice(names, size=min(a["batch_mols"], len(names)), replace=False))
        dd.append(100 * (tr[s]["reward_mean"] - np.mean([et[0][n] for n in b])))
    print(f"**{label}**\n")
    print("- train reward minus start policy on the same molecules (pp), steps 1..: "
          + ", ".join(f"{x:+.1f}" for x in dd))
    vals = []
    for s in sorted(ev):
        d = np.array([ev[s][n] - ev[0][n] for n in ev[0]])
        vals.append(f"{s}: {100 * d.mean():+.2f} ({int((d > 0).sum())}/{len(d)})")
    print("- val (20) paired change vs step 0, pp (# better): " + "; ".join(vals))
    last = max(et)
    if last > 0:
        d = np.array([et[last][n] - et[0][n] for n in et[0]])
        print(f"- train (40) paired change step {last} vs 0: {100 * d.mean():+.2f} pp, better on {int((d > 0).sum())}/40")
    kl = [f"{r['step']}: {r['kl']:.0f}" for s, r in sorted(tr.items()) if (s + 1) % a["ref_update"] == 0]
    print(f"- KL diagnostic (policy vs reference, reference reset every {a['ref_update']} steps), last step of each "
          f"block: " + ", ".join(kl) + "\n")
    return dd


def cost(path, label, n_res):
    """Measured cost per RL step and per evaluation (wall clock from the log's own time stamps)."""
    p = Path(path) / "log.jsonl"
    if not p.exists():
        return
    rows = [json.loads(ln) for ln in open(p)]
    i0 = max(i for i, r in enumerate(rows) if r.get("type") == "eval_val" and r["step"] == 0)
    rows = rows[i0:]
    a = json.load(open(Path(path) / "args.json"))
    ne = a["batch_mols"] * a["group_size"]
    st = [r for r in rows if r.get("type") == "train"]
    ts, tr = np.array([r["t_step"] for r in st]), np.array([r["t_reward"] for r in st])
    print(f"**{label}** ({n_res}; {ne} energies per step)\n")
    print(f"- RL steps: {len(st)}, {ts.sum() / 3600:.2f} h in total; per step median {np.median(ts):.0f} s "
          f"(mean {ts.mean():.0f}); reward phase median {np.median(tr):.0f} s = {np.median(tr) / ne:.1f} s of "
          f"wall clock per energy, {3600 * ne / np.median(tr):.0f} energies per hour")
    ev = []
    for i in range(1, len(rows)):
        r = rows[i]
        if r.get("type") in ("eval_val", "eval_train"):
            ev.append((r["type"], r["step"], r["n"], r["time"] - rows[i - 1]["time"]))
    for kind in ("eval_val", "eval_train"):
        e = [x for x in ev if x[0] == kind]
        if e:
            d = np.array([x[3] for x in e])
            print(f"- {kind} ({e[0][2]} energies, steps {', '.join(str(x[1]) for x in e)}): {d.min():.0f}-{d.max():.0f} s, "
                  f"{d.min() / e[0][2]:.1f}-{d.max() / e[0][2]:.1f} s per energy")
    print(f"- first val eval to last record: {(rows[-1]['time'] - rows[0]['time']) / 3600:.2f} h\n")


def _pids(path):
    p = TD / path
    if not p.exists():
        return set()
    return {int(w) for w in p.read_text().split("\n")[0].replace(":", " ").split() if w.isdigit()}


def gpu_footprint():
    """nvidia-smi footprint of this task's processes on scai7 GPU 4 (gpu_monitor.py samples, every 5 s)."""
    p = TD / "gpu4_monitor_partB.jsonl"
    if not p.exists():
        return
    mon = [json.loads(ln) for ln in open(p) if ln.strip()]
    mon = [m for m in mon if "procs" in m]
    print("| phase | samples | hours | fts MiB median | fts MiB max | max per process | foreign MiB min-max | "
          "GPU util mean |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|")
    for label, pids in (("RL, 4 dm chi-128 workers", _pids("worker_pids.txt")),
                        ("evaluation, 4 v1 chi-256 workers", _pids("v1_worker_pids.txt"))):
        sel = [m for m in mon if pids & {int(k) for k in m["fts"]}]
        if not sel:
            continue
        mine = np.array([sum(v for k, v in m["fts"].items() if int(k) in pids) for m in sel])
        per = max(v for m in sel for k, v in m["fts"].items() if int(k) in pids)
        other = np.array([sum(v for k, v in m["procs"].items() if k not in m["fts"]) for m in sel])
        util = np.array([m["gpu_util"] for m in sel])
        print(f"| {label} | {len(sel)} | {(sel[-1]['t'] - sel[0]['t']) / 3600:.2f} | {np.median(mine):.0f} | "
              f"{mine.max():.0f} | {per} | {other.min():.0f}-{other.max():.0f} | {util.mean():.0f} % |")
    print()


def main():
    print("## RL curves\n")
    curve(RUN, "this run: dm engine chi 128 rewards and evals (scai7 H200, 4 workers x 4 threads)")
    curve(ROOT / "rl_runs/grpo_n29_tn", "Expanse run: v1 zip-up chi 64 rewards and evals (96 CPU workers)")
    print("## drift against the start policy (chi of each run's own reward)\n")
    drift(RUN, "this run (dm chi 128, 6 molecules per step)")
    drift(ROOT / "rl_runs/grpo_n29_tn", "Expanse run (v1 chi 64, 12 molecules per step)")
    print("## cost\n")
    cost(RUN, "this run", "one shared H200, 4 dm chi-128 workers x 4 threads")
    cost(ROOT / "rl_runs/grpo_n29_tn", "Expanse run", "one 128-core node, 96 single-core v1 chi-64 workers")
    print("GPU footprint on scai7 GPU 4 (nvidia-smi, our processes vs the other user's job):\n")
    gpu_footprint()
    val = merge(load(RES / "energy_n29val_base_chi256.json"), load(RES / "energy_n29val_rl4n29_chi256.json"),
                load(NEW / "energy_n29val_rl4n29dm_v1chi256.json"))
    te = merge(load(RES / "energy_n29test_base_chi256.json"), load(RES / "energy_n29test_rl4n29_chi256.json"),
               load(NEW / "energy_n29test_rl4n29dm_v1chi256.json"))
    cands = ["label", "rl4L", "rl4n29", "rl4n29f", "rl4n29dm", "rl4n29dmf"]
    print("## n29 val (20), v1 engine chi 256 (zip margin 1.5), the docs' evaluation setting\n")
    table(val, cands)
    print("## n29 test (24, untouched), v1 engine chi 256\n")
    table(te, cands)
    print("By formula (mean % corr):\n")
    forms = sorted({n.split("_")[0] for n in te})
    print("| formula | n | " + " | ".join(LAB.get(c, c) for c in cands) + " |")
    print("|:--|--:|" + "--:|" * len(cands))
    for fm in forms:
        ns = [n for n in te if n.split("_")[0] == fm]
        print(f"| {fm} | {len(ns)} | " + " | ".join(
            f"{100 * np.mean([te[n][c] for n in ns if c in te[n]]):.1f}" if any(c in te[n] for n in ns) else ""
            for c in cands) + " |")
    print("\n## paired differences (percentage points, mean ± s.e.m., # molecules better)\n")
    out = {}
    for setname, per in (("val20", val), ("test24", te)):
        for a, b in (("rl4n29dm", "rl4n29"), ("rl4n29dmf", "rl4n29f"), ("rl4n29dm", "rl4L"), ("rl4n29dmf", "rl4L"),
                     ("rl4n29dmf", "rl4n29"), ("rl4n29dm", "rl4n29f")):
            pw = pairwise(per, a, b)
            if pw:
                out[f"{setname}:{a}-{b}"] = pw
                print(f"- {setname}: {a} − {b}: {pw['mean_pp']:+.2f} ± {pw['sem_pp']:.2f} pp, better on "
                      f"{pw['better']}/{pw['n']}")
    # the run's own dm chi-128 val evals against v1 chi-256 for the same three policies (step 0 = rl4L)
    rows = [json.loads(ln) for ln in open(RUN / "log.jsonl")] if (RUN / "log.jsonl").exists() else []
    evv = {r["step"]: r["per_mol"] for r in rows if r.get("type") == "eval_val"}
    print()
    for step, tag in ((0, "rl4L"), (20, "rl4n29dm"), (24, "rl4n29dmf")):
        if step in evv and all(tag in val.get(n, {}) for n in evv[step]):
            d = np.array([100 * (evv[step][n] - val[n][tag]) for n in evv[step]])
            print(f"- val20, step {step} ({tag}): dm chi 128 minus v1 chi 256 = {d.mean():+.2f} pp "
                  f"(range {d.min():+.2f} to {d.max():+.2f})")
            out[f"val20:dm128-v1chi256:{tag}"] = {"mean_pp": d.mean(), "min_pp": d.min(), "max_pp": d.max()}
    ctl = load(NEW / "energy_n29val_control_v1chi256.json")
    ref = load(RES / "energy_n29val_rl4n29_chi256.json")
    if ctl and ref:
        d = [(ctl[n][c]["E"] - ref[n][c]["E"]) * 1e3 for n in ctl for c in ctl[n] if c in ref.get(n, {})]
        print(f"\ncontrol (stored v1 chi-256 energies recomputed on this GPU): max |dE| = {max(map(abs, d)):.3f} mHa "
              f"over {len(d)} energies")
        out["control_max_abs_dE_mHa"] = max(map(abs, d))
    json.dump(out, open(NEW / "summary_pairs.json", "w"), indent=1)


if __name__ == "__main__":
    main()
