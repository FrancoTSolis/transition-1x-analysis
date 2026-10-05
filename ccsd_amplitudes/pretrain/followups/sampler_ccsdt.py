#!/usr/bin/env python3
"""CCSD(T) reference for the norb-29 sampler demonstration, computed with pretrain.rl.ccsd_t_refs.one (same RHF,
frozen core and checks as the shared cache) but written to the sampler results only (the shared cache
pretrain/opt_true/results/ccsd_t_refs.json is left untouched).

  OMP_NUM_THREADS=1 python -m pretrain.followups.sampler_ccsdt C3H5N3_rxn4472_P
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pretrain.rl.ccsd_t_refs import one  # noqa: E402

out = ROOT / "pretrain" / "opt_true" / "results" / "followups" / "sampler" / "ccsd_t_refs_sampler.json"
res = json.load(open(out)) if out.exists() else {}
for name in sys.argv[1:]:
    _, r = one(name)
    res[name] = r
    print(name, r, flush=True)
json.dump(res, open(out, "w"), indent=1, sort_keys=True)
