#!/usr/bin/env python3
"""Print a compact progress table for runs_ot/<run>/log.jsonl files.

Usage: python3 -m pretrain.opt_true.summarize runs_ot/fit16_* [--all]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def rows(path):
    try:
        return [json.loads(ln) for ln in open(path) if ln.strip()]
    except FileNotFoundError:
        return []


def main(argv):
    show_all = "--all" in argv
    runs = [a for a in argv if not a.startswith("--")]
    for r in runs:
        p = Path(r) / "log.jsonl"
        rs = rows(p)
        ev = [x for x in rs if "val_resid_median" in x]
        if not ev:
            print(f"{Path(r).name:34s} (no evals yet)")
            continue
        sel = ev if show_all else [ev[0], ev[len(ev) // 2], ev[-1]] if len(ev) > 2 else ev
        done = any(x.get("type") == "done" for x in rs)
        for x in sel:
            print(f"{Path(r).name:34s} ep {x['epoch']:5d}  train resid {x['train_resid_median']:.3f} "
                  f"(imag {x['train_imag_median']:.3f}, obj {x['train_obj_median']:.3f})  "
                  f"val resid {x['val_resid_median']:.3f} (imag {x['val_imag_median']:.3f}, obj {x['val_obj_median']:.3f})"
                  + (f"  hyps {x.get('val_hyp_usage')}" if x.get("val_hyp_usage") else ""))
        if done:
            print(f"{'':34s} done; best val resid {[x for x in rs if x.get('type') == 'done'][0]['best_val_resid_median']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1:])
