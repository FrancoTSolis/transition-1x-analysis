#!/usr/bin/env python
# Copyright 2026 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Restricted UCJ optimization with a second-order BCH energy.

This module implements a real, closed-shell version of the full Newton
algorithm derived in ``ucj_full_newton_updated.tex``.  The objective is the
second-order BCH approximation, not the exact variational UCJ energy.

All many-electron matrix elements are evaluated in the reference-connected
CISD space.  This projection is exact for the expressions used here: every
individual Hamiltonian or UCJ response operator is at most two-body and acts
once on the Hartree--Fock determinant before the expressions are recast as
matrix elements between reference/CISD vectors.
"""

import numpy
import scipy.linalg
import scipy.sparse.linalg

from pyscf import __config__, ao2mo, ci, lib, scf
from pyscf.cc import ccsd
from pyscf.lib import logger


def _transform_tensor4(eri, rotation):
    """Transform a real chemist-notation four-index tensor."""
    return numpy.einsum('ap,bq,cr,ds,abcd->pqrs',
                        rotation, rotation, rotation, rotation, eri,
                        optimize=True)


class _Operator:
    """Real spin-free operator represented by one- and two-electron tensors."""

    def __init__(self, h1, eri):
        self.h1 = numpy.asarray(h1, dtype=float)
        self.eri = numpy.asarray(eri, dtype=float)
        self._eris_cache = None
        self._eref_cache = None

    def copy(self):
        return _Operator(self.h1.copy(), self.eri.copy())

    def scaled(self, factor):
        return _Operator(factor * self.h1, factor * self.eri)

    def add(self, other, alpha=1.0):
        return _Operator(self.h1 + alpha * other.h1,
                         self.eri + alpha * other.eri)

    def norm(self):
        return numpy.sqrt(numpy.vdot(self.h1, self.h1).real
                          + numpy.vdot(self.eri, self.eri).real)

    def rotate(self, rotation):
        h1 = rotation.T @ self.h1 @ rotation
        return _Operator(h1, _transform_tensor4(self.eri, rotation))

    def commute_one_body(self, generator):
        """Return ``[self, K]`` for a real anti-symmetric one-body K."""
        h1 = self.h1 @ generator - generator @ self.h1
        eri = (numpy.einsum('ap,aqrs->pqrs', generator, self.eri,
                            optimize=True)
               + numpy.einsum('bq,pbrs->pqrs', generator, self.eri,
                              optimize=True)
               + numpy.einsum('cr,pqcs->pqrs', generator, self.eri,
                              optimize=True)
               + numpy.einsum('ds,pqrd->pqrs', generator, self.eri,
                              optimize=True))
        return _Operator(h1, eri)


def _jset_operator(g):
    r"""Build ``G = 1/2 sum(pq) g[pq] n[p] n[q]``.

    In the standard spin-free electronic-Hamiltonian representation the
    diagonal one-electron term restores the contraction removed from the
    two-electron excitation operator.
    """
    g = numpy.asarray(g, dtype=float)
    nmo = g.shape[0]
    h1 = numpy.zeros((nmo, nmo))
    h1[numpy.diag_indices(nmo)] = .5 * numpy.diag(g)
    eri = numpy.zeros((nmo,) * 4)
    idx = numpy.arange(nmo)
    eri[idx[:, None], idx[:, None], idx[None, :], idx[None, :]] = g
    return _Operator(h1, eri)


def _q_operator(nmo, p, q):
    """Derivative of G with respect to one packed symmetric g parameter."""
    h1 = numpy.zeros((nmo, nmo))
    eri = numpy.zeros((nmo,) * 4)
    if p == q:
        h1[p, p] = .5
        eri[p, p, p, p] = 1.0
    else:
        eri[p, p, q, q] = 1.0
        eri[q, q, p, p] = 1.0
    return _Operator(h1, eri)


def _operator_eris(op, solver):
    """Adapt an arbitrary Hermitian operator to PySCF's CISD contraction."""
    if op._eris_cache is not None:
        return op._eris_cache, op._eref_cache

    nocc = solver.nocc
    nmo = solver.nmo
    nvir = nmo - nocc
    eri = op.eri
    eris = ccsd._ChemistsERIs(solver.mol)
    eris.nocc = nocc
    eris.mo_coeff = solver.mo_coeff
    eris.fock = op.h1.copy()
    for i in range(nocc):
        eris.fock += 2 * eri[:, :, i, i] - eri[:, i, i, :]

    o = slice(0, nocc)
    v = slice(nocc, nmo)
    eris.oooo = eri[o, o, o, o]
    eris.ovoo = eri[o, v, o, o]
    eris.oovv = eri[o, o, v, v]
    eris.ovvo = eri[o, v, v, o]
    eris.ovov = eri[o, v, o, v]
    ovvv = eri[o, v, v, v].reshape(nocc * nvir, nvir, nvir)
    eris.ovvv = lib.pack_tril(ovvv).reshape(nocc, nvir, -1)
    eris.vvvv = eri[v, v, v, v]

    eref = 2 * numpy.trace(op.h1[o, o])
    for i in range(nocc):
        for j in range(nocc):
            eref += 2 * eri[i, i, j, j] - eri[i, j, j, i]
    op._eris_cache = eris
    op._eref_cache = float(eref)
    return eris, float(eref)


def _pack_g(g, tril_indices=None):
    g = numpy.asarray(g, dtype=float)
    if tril_indices is None:
        tril_indices = numpy.triu_indices(g.shape[0])
    return g[tril_indices].copy()


def _unpack_g(vector, nmo, tril_indices=None):
    if tril_indices is None:
        tril_indices = numpy.triu_indices(nmo)
    vector = numpy.asarray(vector, dtype=float)
    g = numpy.zeros((nmo, nmo))
    g[tril_indices] = vector
    g[(tril_indices[1], tril_indices[0])] = vector
    return g


def _boundary_intersection(step, direction, radius):
    aa = numpy.dot(direction, direction)
    bb = 2 * numpy.dot(step, direction)
    cc = numpy.dot(step, step) - radius ** 2
    disc = max(0.0, bb * bb - 4 * aa * cc)
    return (-bb + numpy.sqrt(disc)) / (2 * aa)


def _steihaug(gradient, hessp, radius, tolerance=None, max_cycle=None):
    """Approximately minimize a quadratic inside a Euclidean trust region."""
    nparam = gradient.size
    if max_cycle is None:
        max_cycle = nparam
    if tolerance is None:
        tolerance = min(.5, numpy.sqrt(numpy.linalg.norm(gradient))) \
                    * numpy.linalg.norm(gradient)
    step = numpy.zeros_like(gradient)
    residual = gradient.copy()
    direction = -residual
    rr = numpy.dot(residual, residual)
    if numpy.sqrt(rr) <= tolerance:
        return step, False, 0

    for cycle in range(1, max_cycle + 1):
        hd = hessp(direction)
        curvature = numpy.dot(direction, hd)
        if curvature <= 0:
            tau = _boundary_intersection(step, direction, radius)
            return step + tau * direction, True, cycle
        alpha = rr / curvature
        trial = step + alpha * direction
        if numpy.linalg.norm(trial) >= radius:
            tau = _boundary_intersection(step, direction, radius)
            return step + tau * direction, True, cycle
        step = trial
        residual = residual + alpha * hd
        rr_new = numpy.dot(residual, residual)
        if numpy.sqrt(rr_new) <= tolerance:
            return step, False, cycle
        beta = rr_new / rr
        direction = -residual + beta * direction
        rr = rr_new
    return step, numpy.linalg.norm(step) >= .99 * radius, max_cycle


class _LocalModel:
    """Energy, gradient, and Hessian action at one accumulated-U origin."""

    def __init__(self, parent, g, rotation):
        self.parent = parent
        self.g = g
        self.rotation = rotation
        self.x = _jset_operator(g).rotate(rotation)
        self.qops = [q.rotate(rotation) for q in parent._qops]
        self.kops = parent._kops
        self.h = parent._hamiltonian
        self.ref = parent._reference
        self._h0 = parent._apply(self.h, self.ref)
        self._x0 = parent._apply(self.x, self.ref)
        self._hx0 = parent._apply(self.h, self._x0)
        self._xx0 = parent._apply(self.x, self._x0)
        self.first_ops = [self.x.commute_one_body(k) for k in self.kops]
        self.first_ops.extend(self.qops)
        self.first_vecs = [parent._apply(op, self.ref)
                           for op in self.first_ops]

        self.orbital_mask = numpy.asarray(
            [op.norm() > parent.lindep for op in self.first_ops[:parent.nkappa]])

    def dot(self, left, right):
        return self.parent._dot(left, right)

    @property
    def energy_elec(self):
        return self.parent.e_ref_elec \
            + self.dot(self._x0, self._hx0) \
            - self.dot(self._h0, self._xx0)

    def gform(self, op, op0=None):
        if op0 is None:
            op0 = self.parent._apply(op, self.ref)
        return (2 * self.dot(op0, self._hx0)
                - self.dot(self._h0, self.parent._apply(op, self._x0))
                - self.dot(self._h0, self.parent._apply(self.x, op0)))

    def bform(self, op_a, vec_a, op_b, vec_b, h_vec_b=None):
        if h_vec_b is None:
            h_vec_b = self.parent._apply(self.h, vec_b)
        return (2 * self.dot(vec_a, h_vec_b)
                - self.dot(self._h0, self.parent._apply(op_a, vec_b))
                - self.dot(self._h0, self.parent._apply(op_b, vec_a)))

    def project(self, vector):
        vector = numpy.asarray(vector, dtype=float).copy()
        vector[:self.parent.nkappa] *= self.orbital_mask
        vg = vector[self.parent.nkappa:]
        phase = self.parent._phase_direction
        vg -= phase * numpy.dot(phase, vg)
        return vector

    def gradient(self):
        values = numpy.asarray([self.gform(op, vec)
                                for op, vec in zip(self.first_ops,
                                                   self.first_vecs)])
        return self.project(values)

    def hessp(self, direction):
        parent = self.parent
        direction = self.project(direction)
        pk = direction[:parent.nkappa]
        pg = direction[parent.nkappa:]

        kdir = parent._kappa_matrix(pk)
        qdir = _jset_operator(_unpack_g(pg, parent.nmo,
                                       parent._g_indices)).rotate(self.rotation)
        adir = self.x.commute_one_body(kdir)
        xdir = adir.add(qdir)
        xdir0 = parent._apply(xdir, self.ref)
        hxdir0 = parent._apply(self.h, xdir0)

        result = numpy.empty(parent.nparam)
        for mu, (op, vec, kmu) in enumerate(zip(
                self.first_ops[:parent.nkappa],
                self.first_vecs[:parent.nkappa], self.kops)):
            value = self.bform(op, vec, xdir, xdir0, hxdir0)
            d1 = op.commute_one_body(kdir)
            d2 = adir.commute_one_body(kmu)
            second = d1.add(d2).scaled(.5)
            second = second.add(qdir.commute_one_body(kmu))
            value += self.gform(second)
            result[mu] = value

        for alpha, (qop, qvec) in enumerate(zip(
                self.qops, self.first_vecs[parent.nkappa:])):
            value = self.bform(qop, qvec, xdir, xdir0, hxdir0)
            second = qop.commute_one_body(kdir)
            value += self.gform(second)
            result[parent.nkappa + alpha] = value
        return self.project(result)

    def kappa_hessp(self, direction):
        """Orbital-orbital Hessian action used to construct the initial seed."""
        parent = self.parent
        direction = numpy.asarray(direction, dtype=float).copy()
        direction *= self.orbital_mask
        kdir = parent._kappa_matrix(direction)
        adir = self.x.commute_one_body(kdir)
        adir0 = parent._apply(adir, self.ref)
        hadir0 = parent._apply(self.h, adir0)

        result = numpy.empty(parent.nkappa)
        for mu, (op, vec, kmu) in enumerate(zip(
                self.first_ops[:parent.nkappa],
                self.first_vecs[:parent.nkappa], self.kops)):
            value = self.bform(op, vec, adir, adir0, hadir0)
            d1 = op.commute_one_body(kdir)
            d2 = adir.commute_one_body(kmu)
            value += self.gform(d1.add(d2).scaled(.5))
            result[mu] = value
        result *= self.orbital_mask
        return result

    def g_hessp(self, direction):
        """Pure Jastrow Hessian action, projected out of the phase gauge."""
        parent = self.parent
        direction = numpy.asarray(direction, dtype=float).copy()
        phase = parent._phase_direction
        direction -= phase * numpy.dot(phase, direction)
        qdir = _jset_operator(_unpack_g(
            direction, parent.nmo, parent._g_indices)).rotate(self.rotation)
        qdir0 = parent._apply(qdir, self.ref)
        hqdir0 = parent._apply(self.h, qdir0)

        result = numpy.asarray([
            self.bform(qop, qvec, qdir, qdir0, hqdir0)
            for qop, qvec in zip(
                self.qops, self.first_vecs[parent.nkappa:])])
        result -= phase * numpy.dot(phase, result)
        return result

    def _mixed_kappa_g(self, direction):
        """Return B times a packed Jastrow direction."""
        full = numpy.zeros(self.parent.nparam)
        full[self.parent.nkappa:] = direction
        return self.hessp(full)[:self.parent.nkappa]

    def _mixed_g_kappa(self, direction):
        """Return B.T times an occupied--virtual direction."""
        full = numpy.zeros(self.parent.nparam)
        full[:self.parent.nkappa] = direction
        return self.hessp(full)[self.parent.nkappa:]

    def _solve_shifted_kappa(self, rhs, level_shift, tolerance, max_cycle):
        """Apply (H_kk + level_shift I)^-1 in the active kappa space."""
        rhs = numpy.asarray(rhs, dtype=float)
        indices = numpy.flatnonzero(self.orbital_mask)
        if indices.size == 0 or numpy.linalg.norm(rhs[indices]) == 0:
            return numpy.zeros(self.parent.nkappa)

        def action(vector):
            full = numpy.zeros(self.parent.nkappa)
            full[indices] = vector
            return (self.kappa_hessp(full)[indices]
                    + level_shift * vector)

        operator = scipy.sparse.linalg.LinearOperator(
            (indices.size, indices.size), matvec=action, dtype=float)
        solution, info = scipy.sparse.linalg.minres(
            operator, rhs[indices], rtol=tolerance, maxiter=max_cycle)
        if info != 0:
            logger.warn(self.parent, 'Shifted H_kk MINRES returned info = %s',
                        info)
        full = numpy.zeros(self.parent.nkappa)
        full[indices] = solution
        return full

    def relaxed_g_hop(self, kappa_shift=1e-3, g_shift=0.0,
                      tolerance=1e-7, max_cycle=None):
        r"""Return the orbital-relaxed g gradient and Schur action.

        The returned matrix-free action is

        ``C + g_shift*I - B.T (A + kappa_shift*I)^-1 B``.
        """
        parent = self.parent
        if max_cycle is None:
            max_cycle = max(20, 2 * parent.nkappa)
        gradient = self.gradient()
        rk = gradient[:parent.nkappa]
        rg = gradient[parent.nkappa:]

        # At X=0, A, B, and r_k vanish identically.  Avoid an unnecessary
        # singular orbital solve; the relaxed block then reduces exactly to C.
        no_orbital_response = self.x.norm() <= parent.lindep
        if no_orbital_response:
            relaxed_gradient = rg.copy()
        else:
            response = self._solve_shifted_kappa(
                rk, kappa_shift, tolerance, max_cycle)
            relaxed_gradient = rg - self._mixed_g_kappa(response)

        phase = parent._phase_direction
        relaxed_gradient -= phase * numpy.dot(phase, relaxed_gradient)

        def schur(direction):
            direction = numpy.asarray(direction, dtype=float).copy()
            direction -= phase * numpy.dot(phase, direction)
            result = self.g_hessp(direction) + g_shift * direction
            if not no_orbital_response:
                mixed = self._mixed_kappa_g(direction)
                response = self._solve_shifted_kappa(
                    mixed, kappa_shift, tolerance, max_cycle)
                result -= self._mixed_g_kappa(response)
            result -= phase * numpy.dot(phase, result)
            return result

        return relaxed_gradient, schur


class RUCJ(lib.StreamObject):
    """Real closed-shell UCJ-BCH2 optimizer.

    Args:
        mf : converged :class:`pyscf.scf.hf.RHF` object

    The reported ``e_tot`` is a second-order BCH energy and is not a strict
    variational upper bound.
    """

    max_cycle = getattr(__config__, 'ucj_ucj_RUCJ_max_cycle', 50)
    max_reject = getattr(__config__, 'ucj_ucj_RUCJ_max_reject', 10)
    conv_tol_grad = getattr(__config__, 'ucj_ucj_RUCJ_conv_tol_grad', 1e-6)
    conv_tol_step = getattr(__config__, 'ucj_ucj_RUCJ_conv_tol_step', 1e-6)
    conv_tol_energy = getattr(__config__, 'ucj_ucj_RUCJ_conv_tol_energy', 1e-9)
    trust_radius = getattr(__config__, 'ucj_ucj_RUCJ_trust_radius', .25)
    trust_radius_max = getattr(__config__, 'ucj_ucj_RUCJ_trust_radius_max', 1.0)
    trust_accept = getattr(__config__, 'ucj_ucj_RUCJ_trust_accept', .1)
    lindep = getattr(__config__, 'ucj_ucj_RUCJ_lindep', 1e-10)
    kappa_seed = getattr(__config__, 'ucj_ucj_RUCJ_kappa_seed', 1e-2)
    g_seed_radius = getattr(__config__, 'ucj_ucj_RUCJ_g_seed_radius', 5e-2)
    relaxed_kappa_shift = getattr(
        __config__, 'ucj_ucj_RUCJ_relaxed_kappa_shift', 1e-3)
    relaxed_g_shift = getattr(
        __config__, 'ucj_ucj_RUCJ_relaxed_g_shift', 0.0)
    seed_eig_tol = getattr(__config__, 'ucj_ucj_RUCJ_seed_eig_tol', 1e-5)
    seed_eig_max_cycle = getattr(
        __config__, 'ucj_ucj_RUCJ_seed_eig_max_cycle', 40)

    _keys = {
        'max_cycle', 'max_reject', 'conv_tol_grad', 'conv_tol_step',
        'conv_tol_energy', 'trust_radius', 'trust_radius_max',
        'trust_accept', 'lindep', 'converged', 'e_bch2', 'e_tot',
        'kappa_seed', 'g_seed_radius', 'relaxed_kappa_shift',
        'relaxed_g_shift', 'seed_eig_tol', 'seed_eig_max_cycle',
        'e_corr', 'g', 'orbital_rotation', 'history',
        'seed_eigenvalue', 'seed_energy_change', 'seed_sign',
        'kappa_seed_vector', 'g_seed_eigenvalue', 'g_seed_vector',
        'g_seed_energy_change', 'cisd_e_corr', 'cisd_c0',
        'cisd_singles_norm',
    }

    def __init__(self, mf):
        if isinstance(mf, (scf.uhf.UHF, scf.rohf.ROHF, scf.ghf.GHF)) \
                or not isinstance(mf, scf.hf.RHF):
            raise TypeError('RUCJ requires a restricted closed-shell RHF object')
        if not mf.converged:
            raise RuntimeError('The RHF reference must be converged')
        if numpy.iscomplexobj(mf.mo_coeff):
            raise NotImplementedError('Complex orbitals are not supported')
        mo_occ = numpy.asarray(mf.mo_occ)
        if not numpy.all(numpy.logical_or(mo_occ == 0, mo_occ == 2)):
            raise ValueError('RUCJ requires closed-shell occupations of 2 or 0')

        self._scf = mf
        self.mol = mf.mol
        self.stdout = mf.stdout
        self.verbose = mf.verbose
        self.max_memory = mf.max_memory
        self.mo_coeff = numpy.asarray(mf.mo_coeff, dtype=float)
        self.mo_occ = mo_occ.copy()
        self.nmo = self.mo_coeff.shape[1]
        self.nocc = int(numpy.count_nonzero(self.mo_occ > 0))
        self.nvir = self.nmo - self.nocc
        self._g_indices = numpy.triu_indices(self.nmo)
        self.ng = len(self._g_indices[0])
        self.nkappa = self.nocc * self.nvir
        self.nparam = self.nkappa + self.ng

        self._kops = []
        for a in range(self.nocc, self.nmo):
            for i in range(self.nocc):
                k = numpy.zeros((self.nmo, self.nmo))
                k[a, i] = 1.0
                k[i, a] = -1.0
                self._kops.append(k)
        self._qops = [_q_operator(self.nmo, int(p), int(q))
                      for p, q in zip(*self._g_indices)]
        phase = numpy.ones(self.ng)
        self._phase_direction = phase / numpy.linalg.norm(phase)

        hcore = mf.get_hcore()
        h1 = self.mo_coeff.T @ hcore @ self.mo_coeff
        eri = ao2mo.restore(1, ao2mo.kernel(mf.mol, self.mo_coeff), self.nmo)
        self._hamiltonian = _Operator(h1, eri)
        self._ci = ci.CISD(mf)
        self._ci.nmo = self.nmo
        self._ci.nocc = self.nocc
        self._ci.direct = False
        self._reference = numpy.zeros(self._ci.vector_size())
        self._reference[0] = 1.0
        self.e_ref_elec = self._dot(
            self._reference, self._apply(self._hamiltonian, self._reference))
        self.e_nuc = self.mol.energy_nuc()

        self.converged = False
        self.e_bch2 = None
        self.e_tot = None
        self.e_corr = None
        self.g = None
        self.orbital_rotation = None
        self.history = []
        self._last_model = None
        self.seed_eigenvalue = None
        self.seed_energy_change = None
        self.seed_sign = None
        self.kappa_seed_vector = None
        self.g_seed_eigenvalue = None
        self.g_seed_vector = None
        self.g_seed_energy_change = None
        self.cisd_e_corr = None
        self.cisd_c0 = None
        self.cisd_singles_norm = None

    def _apply(self, op, vector):
        eris, eref = _operator_eris(op, self)
        return self._ci.contract(vector, eris) + eref * vector

    def _dot(self, left, right):
        return float(self._ci._dot(left, right, self.nmo, self.nocc))

    def _kappa_matrix(self, vector):
        matrix = numpy.zeros((self.nmo, self.nmo))
        if self.nkappa:
            block = numpy.asarray(vector).reshape(self.nvir, self.nocc)
            matrix[self.nocc:, :self.nocc] = block
            matrix[:self.nocc, self.nocc:] = -block.T
        return matrix

    def dump_flags(self, verbose=None):
        log = logger.new_logger(self, verbose)
        log.info('')
        log.info('******** RHF UCJ-BCH2 Newton solver ********')
        log.info('nocc = %d  nvir = %d  nkappa = %d  ng = %d',
                 self.nocc, self.nvir, self.nkappa, self.ng)
        log.info('max_cycle = %d  max_reject = %d',
                 self.max_cycle, self.max_reject)
        log.info('conv_tol_grad = %g  conv_tol_step = %g',
                 self.conv_tol_grad, self.conv_tol_step)
        log.info('conv_tol_energy = %g', self.conv_tol_energy)
        log.info('trust radius = %g  max = %g',
                 self.trust_radius, self.trust_radius_max)
        log.info('CISD kappa seed norm = %g  relaxed-Hgg g radius = %g',
                 self.kappa_seed, self.g_seed_radius)
        log.info('relaxed-Hgg shifts: kappa = %g  g = %g',
                 self.relaxed_kappa_shift, self.relaxed_g_shift)
        log.warn('UCJ-BCH2 is a truncated BCH energy, not a variational bound')
        return self

    def _seed_orbital_mask(self):
        """Return symmetry-allowed occupied-virtual variables when available."""
        mask = numpy.ones(self.nkappa, dtype=bool)
        if not self.mol.symmetry:
            return mask
        try:
            from pyscf import symm
            orbsym = symm.label_orb_symm(
                self.mol, self.mol.irrep_id, self.mol.symm_orb,
                self.mo_coeff)
        except (AttributeError, ValueError):
            logger.warn(self, 'Unable to label orbital symmetries; seeding all '
                        'occupied-virtual rotations')
            return mask
        allowed = []
        for a in range(self.nocc, self.nmo):
            for i in range(self.nocc):
                allowed.append(orbsym[a] == orbsym[i])
        return numpy.asarray(allowed, dtype=bool)

    def _energy_from_model(self, model):
        return model.energy_elec + self.e_nuc

    def _initial_rotation(self):
        """Build a small occupied--virtual rotation from CISD singles."""
        identity = numpy.eye(self.nmo)
        if self.nkappa == 0 or self.kappa_seed <= 0:
            return identity

        seed_ci = ci.CISD(self._scf)
        seed_ci.verbose = max(0, self.verbose - 1)
        e_corr, civec = seed_ci.kernel()
        c0, c1, _ = seed_ci.cisdvec_to_amplitudes(civec)
        self.cisd_e_corr = float(e_corr)
        self.cisd_c0 = float(c0)
        if abs(c0) <= self.lindep:
            logger.warn(self, 'CISD reference coefficient is too small for a '
                        'singles orbital seed')
            return identity

        # PySCF stores c1 as (i,a), whereas kappa is packed with a outermost
        # and i innermost.  c1/c0 is the first-order Thouless direction.
        direction = numpy.asarray(c1 / c0, dtype=float).T.ravel()
        direction *= self._seed_orbital_mask()
        self.cisd_singles_norm = float(numpy.linalg.norm(direction))
        if self.cisd_singles_norm <= self.lindep:
            logger.warn(self, 'CISD singles are zero; no orbital seed is '
                        'available for the relaxed-Hgg initialization')
            return identity
        direction /= self.cisd_singles_norm
        kappa = self.kappa_seed * direction
        self.kappa_seed_vector = kappa.copy()
        self.seed_sign = 1
        self.seed_eigenvalue = None
        self.seed_energy_change = 0.0
        logger.info(self, 'Initial CISD singles seed: |c1/c0|2 = %.6g  '
                    '|kappa|2 = %.3g', self.cisd_singles_norm,
                    self.kappa_seed)
        return scipy.linalg.expm(self._kappa_matrix(kappa))

    def _phase_complement(self):
        """Orthonormal basis for packed-g space excluding the N^2 gauge."""
        return scipy.linalg.null_space(self._phase_direction[numpy.newaxis])

    def get_init_guess(self, rotation=None):
        """Predict initial g from the lowest orbital-relaxed H_gg mode.

        A finite CISD-singles orbital seed is used when ``rotation`` is not
        supplied.  At ``g=0`` the mixed orbital block vanishes structurally,
        so the Schur complement correctly reduces to the bare ``H_gg`` at
        that finite orientation.  A negative mode is followed to the bounded
        seed radius; if no negative mode exists, the quadratic prediction is
        the zero vector.
        """
        if rotation is None:
            rotation = self._initial_rotation()
        else:
            rotation = self._validate_u(rotation)
        zero_g = numpy.zeros((self.nmo, self.nmo))
        model = _LocalModel(self, zero_g, rotation)
        relaxed_gradient, relaxed_hessp = model.relaxed_g_hop(
            self.relaxed_kappa_shift, self.relaxed_g_shift)
        basis = self._phase_complement()
        if basis.shape[1] == 0 or self.g_seed_radius <= 0:
            self.g_seed_eigenvalue = 0.0
            self.g_seed_vector = numpy.zeros(self.ng)
            self.g_seed_energy_change = 0.0
            return zero_g

        def reduced_hessp(vector):
            return basis.T @ relaxed_hessp(basis @ vector)

        dimension = basis.shape[1]
        if dimension == 1:
            eigenvector = numpy.ones(1)
            eigenvalue = float(reduced_hessp(eigenvector)[0])
        else:
            operator = scipy.sparse.linalg.LinearOperator(
                (dimension, dimension), matvec=reduced_hessp, dtype=float)
            v0 = numpy.linspace(1.0, 2.0, dimension)
            v0 /= numpy.linalg.norm(v0)
            try:
                values, vectors = scipy.sparse.linalg.eigsh(
                    operator, k=1, which='SA', v0=v0,
                    tol=self.seed_eig_tol,
                    maxiter=self.seed_eig_max_cycle,
                    ncv=min(dimension, 12))
                eigenvalue = float(values[0])
                eigenvector = vectors[:, 0]
            except scipy.sparse.linalg.ArpackNoConvergence as err:
                if err.eigenvalues is None or len(err.eigenvalues) == 0:
                    logger.warn(self, 'Relaxed-Hgg eigenproblem did not '
                                'converge; using its deterministic trial mode')
                    eigenvector = v0
                    eigenvalue = float(numpy.dot(
                        eigenvector, reduced_hessp(eigenvector)))
                else:
                    eigenvalue = float(err.eigenvalues[0])
                    eigenvector = err.eigenvectors[:, 0]
                    logger.warn(self, 'Using a partially converged relaxed-Hgg '
                                'eigenpair')

        direction = basis @ eigenvector
        direction /= numpy.linalg.norm(direction)
        pivot = numpy.argmax(numpy.abs(direction))
        if direction[pivot] < 0:
            direction *= -1
        if (eigenvalue < -self.lindep
                or numpy.linalg.norm(relaxed_gradient) > self.lindep):
            # With zero relaxed gradient, an indefinite quadratic is minimized
            # on the trust boundary along its lowest-curvature eigenvector.
            packed = self.g_seed_radius * direction
        else:
            packed = numpy.zeros(self.ng)
        packed -= self._phase_direction * numpy.dot(
            self._phase_direction, packed)
        g = _unpack_g(packed, self.nmo, self._g_indices)
        self.g_seed_eigenvalue = eigenvalue
        self.g_seed_vector = packed.copy()
        self.g_seed_energy_change = (
            self._energy_from_model(_LocalModel(self, g, rotation))
            - self._energy_from_model(model))
        logger.info(self, 'Initial relaxed-Hgg mode: curvature = %.6g  '
                    '|g|2 = %.3g  dE = %.6g', eigenvalue,
                    numpy.linalg.norm(packed), self.g_seed_energy_change)
        return g

    def _validate_g(self, g):
        g = numpy.asarray(g, dtype=float)
        if g.shape == (self.ng,):
            g = _unpack_g(g, self.nmo, self._g_indices)
        if g.shape != (self.nmo, self.nmo):
            raise ValueError('g must be a symmetric (nmo,nmo) matrix or packed vector')
        if not numpy.allclose(g, g.T, atol=1e-12):
            raise ValueError('g must be symmetric')
        return .5 * (g + g.T)

    def _validate_u(self, rotation):
        rotation = numpy.asarray(rotation, dtype=float)
        if rotation.shape != (self.nmo, self.nmo):
            raise ValueError('u0 must have shape (nmo,nmo)')
        if not numpy.allclose(rotation.T @ rotation, numpy.eye(self.nmo),
                             atol=1e-10):
            raise ValueError('u0 must be orthogonal')
        return rotation.copy()

    def gen_g_hop(self, g=None, u=None):
        """Return the projected gradient and matrix-free Hessian action."""
        if g is None:
            if self.g is None:
                raise RuntimeError('No Jastrow parameters are available')
            g = self.g
        if u is None:
            u = self.orbital_rotation
            if u is None:
                u = numpy.eye(self.nmo)
        model = _LocalModel(self, self._validate_g(g), self._validate_u(u))
        self._last_model = model
        return model.gradient(), model.hessp

    def gen_relaxed_g_hop(self, g=None, u=None, kappa_shift=None,
                          g_shift=None):
        """Return the relaxed g gradient and Schur-complement action."""
        if g is None:
            if self.g is None:
                raise RuntimeError('No Jastrow parameters are available')
            g = self.g
        if u is None:
            u = self.orbital_rotation
            if u is None:
                u = numpy.eye(self.nmo)
        if kappa_shift is None:
            kappa_shift = self.relaxed_kappa_shift
        if g_shift is None:
            g_shift = self.relaxed_g_shift
        model = _LocalModel(self, self._validate_g(g), self._validate_u(u))
        self._last_model = model
        return model.relaxed_g_hop(kappa_shift, g_shift)

    def energy(self, g=None, u=None):
        """Evaluate the total UCJ-BCH2 energy for ``g`` and ``u``."""
        if g is None:
            g = self.g
        if u is None:
            u = self.orbital_rotation
        if g is None or u is None:
            raise RuntimeError('g and orbital rotation have not been initialized')
        model = _LocalModel(self, self._validate_g(g), self._validate_u(u))
        return model.energy_elec + self.e_nuc

    def kernel(self, g0=None, u0=None):
        log = logger.new_logger(self)
        self.dump_flags()
        if u0 is None:
            rotation = self._initial_rotation()
        else:
            rotation = self._validate_u(u0)
            self.seed_eigenvalue = None
            self.seed_energy_change = None
            self.seed_sign = None
            self.kappa_seed_vector = None
        if g0 is None:
            g = self.get_init_guess(rotation)
        else:
            g = self._validate_g(g0)
            self.g_seed_eigenvalue = None
            self.g_seed_vector = None
            self.g_seed_energy_change = None
        radius = min(float(self.trust_radius), float(self.trust_radius_max))
        self.history = []
        self.converged = False
        last_step = numpy.inf
        last_de = numpy.inf

        for macro in range(1, self.max_cycle + 1):
            model = _LocalModel(self, g, rotation)
            self._last_model = model
            energy = model.energy_elec + self.e_nuc
            gradient = model.gradient()
            norm_grad = numpy.linalg.norm(gradient, numpy.inf)

            # A stationary starting point has the zero Newton step and zero
            # energy change by definition.  Avoid sending that case through
            # the trust solver merely because no previous step exists.
            if norm_grad < self.conv_tol_grad and not self.history:
                last_step = 0.0
                last_de = 0.0
                self.converged = True
                break
            if (norm_grad < self.conv_tol_grad
                    and last_step < self.conv_tol_step
                    and abs(last_de) < self.conv_tol_energy):
                self.converged = True
                break

            accepted = False
            for attempt in range(1, self.max_reject + 1):
                step, at_boundary, inner = _steihaug(
                    gradient, model.hessp, radius,
                    max_cycle=min(self.nparam, 50))
                hstep = model.hessp(step)
                predicted = -(numpy.dot(gradient, step)
                              + .5 * numpy.dot(step, hstep))
                if predicted <= 0 or not numpy.isfinite(predicted):
                    radius *= .25
                    continue

                pk = step[:self.nkappa]
                pg = step[self.nkappa:]
                trial_g = g + _unpack_g(pg, self.nmo, self._g_indices)
                trial_u = rotation @ scipy.linalg.expm(self._kappa_matrix(pk))
                trial_model = _LocalModel(self, trial_g, trial_u)
                trial_energy = trial_model.energy_elec + self.e_nuc
                actual = energy - trial_energy
                ratio = actual / predicted

                if ratio < .25:
                    radius *= .25
                elif ratio > .75 and at_boundary:
                    radius = min(2 * radius, self.trust_radius_max)

                if actual > 0 and ratio >= self.trust_accept:
                    g = trial_g
                    rotation = trial_u
                    accepted = True
                    last_step = numpy.linalg.norm(step, numpy.inf)
                    kappa_step_norm = numpy.linalg.norm(pk, numpy.inf)
                    g_step_norm = numpy.linalg.norm(pg, numpy.inf)
                    last_de = trial_energy - energy
                    self.history.append({
                        'cycle': macro, 'energy': trial_energy,
                        'gradient_norm': norm_grad,
                        'step_norm': last_step, 'energy_change': last_de,
                        'kappa_step_norm': kappa_step_norm,
                        'g_step_norm': g_step_norm,
                        'trust_radius': radius, 'ratio': ratio,
                        'inner_cycles': inner, 'rejections': attempt - 1,
                    })
                    log.info('UCJ cycle %d  E(BCH2) = %.15g  dE = %.3g  '
                             '|grad|inf = %.3g  |pk|inf = %.3g  '
                             '|pg|inf = %.3g  rho = %.3g',
                             macro, trial_energy, last_de, norm_grad,
                             kappa_step_norm, g_step_norm, ratio)
                    break

            if not accepted:
                log.warn('UCJ trust-region step failed after %d attempts',
                         self.max_reject)
                break

        final_model = _LocalModel(self, g, rotation)
        final_energy = final_model.energy_elec + self.e_nuc
        final_gradient = final_model.gradient()
        if (numpy.linalg.norm(final_gradient, numpy.inf) < self.conv_tol_grad
                and last_step < self.conv_tol_step
                and abs(last_de) < self.conv_tol_energy):
            self.converged = True

        self.g = g
        self.orbital_rotation = rotation
        self.e_bch2 = self.e_tot = final_energy
        self.e_corr = final_energy - self._scf.e_tot
        self._last_model = final_model
        if self.converged:
            log.note('UCJ-BCH2 converged  E = %.15g  E_corr = %.15g',
                     self.e_tot, self.e_corr)
        else:
            log.note('UCJ-BCH2 not converged  E = %.15g  E_corr = %.15g',
                     self.e_tot, self.e_corr)
        return self.e_corr, self.g, self.orbital_rotation


def UCJ(mf):
    """Create a restricted UCJ-BCH2 solver from a converged RHF object."""
    return RUCJ(mf)


__all__ = ['UCJ', 'RUCJ']
