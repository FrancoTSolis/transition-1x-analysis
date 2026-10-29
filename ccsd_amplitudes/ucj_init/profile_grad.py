import pickle, time, sys, numpy as np, torch
from ucj_init.lucj_grad_gpu import LUCJGradGPU, x_from_uz
name = sys.argv[1]
d = pickle.load(open("runs_ot/energy_tasks/task1_sqd_vs_lin/cands_smallval.pkl", "rb"))
U, Z, t1 = [(U, Z, t1) for k, n, U, Z, t1 in d["tasks"] if n == name and k[1] == "label"][0]
ge = LUCJGradGPU.from_npz(name, t1=t1, verbose=True)
x = x_from_uz(np.asarray(U), np.asarray(Z))
ge.energy_grad(x)
sync = torch.cuda.synchronize
def tm(f, *a, n=3):
    sync(); t0 = time.time()
    for _ in range(n): f(*a)
    sync(); return (time.time() - t0) / n
phi, lam, sig = ge.phi, ge.lam, ge.sig
W = np.linalg.qr(np.random.default_rng(0).normal(size=(ge.norb, ge.norb)) + 1j*np.random.default_rng(1).normal(size=(ge.norb, ge.norb)))[0]
print("rotate", tm(ge._rotate, phi, W))
print("phase", tm(ge._phase, phi, np.zeros(ge.norb-1), np.zeros(ge.norb)))
print("sigma", tm(ge.sigma, phi, sig))
print("trdm1", tm(ge._trdm1, lam, phi))
print("phase_grads", tm(ge._phase_grads, lam, phi))
print("total", tm(ge.energy_grad, x))
print("mem GB", torch.cuda.max_memory_allocated()/1e9)
