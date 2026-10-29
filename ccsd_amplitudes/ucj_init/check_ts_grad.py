"""Check LUCJGradGPU on a real TS molecule against the verified energy engine (pretrain.rl.gpu_energy):
energy of the stored candidates, and directional derivatives by central differences of engine energies."""
import pickle, sys, time
import numpy as np
from pretrain.rl.gpu_energy import LUCJEnergyGPU
from ucj_init.lucj_grad_gpu import LUCJGradGPU, x_from_uz, uz_from_x

name = sys.argv[1] if len(sys.argv) > 1 else "C2H3N_rxn2858_TS"
d = pickle.load(open("runs_ot/energy_tasks/task1_sqd_vs_lin/cands_smallval.pkl", "rb"))
cands = {k[1]: (U, Z, t1) for k, n, U, Z, t1 in d["tasks"] if n == name}
t1 = next(iter(cands.values()))[2]
t0 = time.time()
ge = LUCJGradGPU.from_npz(name, t1=t1, verbose=True)
print(f"grad engine build {time.time()-t0:.1f}s")
ref = LUCJEnergyGPU.from_npz(name, max_mem_gb=3.0)
for tag, (U, Z, tt) in cands.items():
    x = x_from_uz(np.asarray(U), np.asarray(Z))
    U2, Z2 = uz_from_x(x, ge.norb, 2)
    E = ge.energy_grad(x, want_grad=False)
    Er = ref.energy(U2, Z2, t1=tt)
    Eo = ref.energy(np.asarray(U), np.asarray(Z), t1=tt)
    print(f"{tag:9s} grad-engine {E:.8f}  engine(x) {Er:.8f}  diff {E-Er:+.1e}  engine(orig U,Z) {Eo:.8f}  corr% {ge.corr_pct(E):6.1f}")
U, Z, tt = cands["label"]
x = x_from_uz(np.asarray(U), np.asarray(Z))
t0 = time.time(); E, g = ge.energy_grad(x); tg = time.time() - t0
t0 = time.time(); E, g = ge.energy_grad(x); tg = time.time() - t0
print(f"energy+grad {tg:.2f}s  |g| {np.linalg.norm(g):.3e}")
rng = np.random.default_rng(0)
for h in (1e-3,):
    for _ in range(4):
        dvec = rng.normal(size=len(x)); dvec /= np.linalg.norm(dvec)
        Up, Zp = uz_from_x(x + h * dvec, ge.norb, 2); Um, Zm = uz_from_x(x - h * dvec, ge.norb, 2)
        fd = (ref.energy(Up, Zp, t1=tt) - ref.energy(Um, Zm, t1=tt)) / (2 * h)
        print(f"  dir deriv: grad {g @ dvec:+.6e}  fd(engine) {fd:+.6e}")
    dvec = g / np.linalg.norm(g)
    Up, Zp = uz_from_x(x + h * dvec, ge.norb, 2); Um, Zm = uz_from_x(x - h * dvec, ge.norb, 2)
    fd = (ref.energy(Up, Zp, t1=tt) - ref.energy(Um, Zm, t1=tt)) / (2 * h)
    print(f"  along g: grad {g @ dvec:+.6e}  fd(engine) {fd:+.6e}")
