"""One TN (MPS) LUCJ energy for an n29 molecule from the RL policy's output (single-threaded), for cost on Expanse."""
import argparse, json, os, sys, time
import numpy as np, torch
sys.path.insert(0, os.getcwd())
ap = argparse.ArgumentParser(); ap.add_argument("--name"); ap.add_argument("--chi", type=int, default=64); ap.add_argument("--out")
a = ap.parse_args()
torch.set_num_threads(1)
from pretrain.opt_true.eval_slot import load_slot
from pretrain.rl.grpo_slot import Mol, Policy, to_flat
from pretrain.rl.tn_energy import LUCJEnergyTN
import copy
idx = json.load(open("rhf_dataset/_index.json"))
slot, a_, _, _ = load_slot("runs_ot/slotall_T4/best.pt", "cpu")
frozen = copy.deepcopy(slot).eval()
pol = Policy(slot, a_["d"]); ck = torch.load("rl_runs/grpo_slot_T4p/policy_best.pt", map_location="cpu", weights_only=False); pol.load_state_dict(ck["policy"]); pol.eval()
t0 = time.time()
m = Mol(a.name, idx, "cpu", 0.005, frozen, 3)
with torch.no_grad():
    dK, dZ = pol(m.x, 0.005); U, Z, r = m.realize(to_flat(dK.double(), dZ.double(), m.fi))
t1 = time.time()
h = np.load(f"rhf_hamiltonians/{a.name}.npz")
scratch = os.path.join(os.environ.get("TMPDIR", "/tmp"), f"tn_{a.name}_{os.getpid()}"); os.makedirs(scratch, exist_ok=True)
ev = LUCJEnergyTN(h["one_body"], h["two_body"], float(h["constant"]), int(h["norb"]), (int(h["nelec_a"]), int(h["nelec_b"])),
                  max_bond=a.chi, device="cpu", name=a.name, block2_threads=1, scratch=scratch)
t2 = time.time()
E, info = ev.energy(U, Z, t1=m.t1)
t3 = time.time()
cf = (float(h["e_hf"]) - E) / (float(h["e_hf"]) - float(h["e_ccsd"]))
rec = {"name": a.name, "chi": a.chi, "E": E, "corr_pct": 100 * cf, "resid": r, "t_policy": t1 - t0, "t_setup": t2 - t1, "t_energy": t3 - t2,
       "info": {k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v)) for k, v in info.items()}}
print(json.dumps(rec), flush=True)
if a.out: open(a.out, "a").write(json.dumps(rec) + "\n")
