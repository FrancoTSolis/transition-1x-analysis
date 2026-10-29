#!/usr/bin/env python
"""Compare RHF, RMP2, RCCSD, and UCJ-BCH2 for LiH/6-31G.

UCJ-BCH2 is a truncated BCH objective and is not a variational upper bound;
the energies below should therefore be compared as method results, not by an
assumed ordering.
"""

import numpy

from pyscf import cc, gto, mp, scf, ucj


mol = gto.M(
    atom='Li 0 0 0; H 0 0 1.60',
    basis='6-31g',
    unit='Angstrom',
    charge=0,
    spin=0,
    symmetry=False,
    verbose=4,
)

mf = scf.RHF(mol).run(conv_tol=1e-10)
pt = mp.MP2(mf).run()
mycc = cc.CCSD(mf).run(conv_tol=1e-9, max_cycle=100)
myucj = ucj.UCJ(mf).run()

rows = [
    ('RHF', mf.e_tot, 0.0, mf.converged),
    ('RMP2', pt.e_tot, pt.e_corr, True),
    ('RCCSD', mycc.e_tot, mycc.e_corr, mycc.converged),
    ('UCJ-BCH2', myucj.e_tot, myucj.e_corr, myucj.converged),
]

print('\nLiH at 1.60 Angstrom, 6-31G, all electrons')
print(f"{'Method':12s} {'Total energy / Eh':>22s} "
      f"{'Correlation / Eh':>22s} {'Converged':>11s}")
for name, total, correlation, converged in rows:
    print(f'{name:12s} {total:22.14f} {correlation:22.14f} '
          f'{str(bool(converged)):>11s}')

final_gradient, _ = myucj.gen_g_hop()
print('\nUCJ macroiterations:', len(myucj.history))
print('UCJ CISD singles norm:', myucj.cisd_singles_norm)
print('UCJ CISD kappa seed norm:', myucj.kappa_seed)
print('UCJ relaxed-Hgg seed curvature:', myucj.g_seed_eigenvalue)
print('UCJ relaxed-Hgg g seed norm:',
      None if myucj.g_seed_vector is None else
      numpy.linalg.norm(myucj.g_seed_vector))
print('UCJ relaxed-Hgg seed energy change:', myucj.g_seed_energy_change)
print('UCJ final projected |gradient|_inf:',
      numpy.linalg.norm(final_gradient, numpy.inf))
if myucj.history:
    print('UCJ final |step|_inf:', myucj.history[-1]['step_norm'])
    print('UCJ final |kappa step|_inf:',
          myucj.history[-1]['kappa_step_norm'])
    print('UCJ final |g step|_inf:', myucj.history[-1]['g_step_norm'])
    print('UCJ final |energy change|:',
          abs(myucj.history[-1]['energy_change']))
