"""Unitary cluster--Jastrow methods.

The implementation currently provides the real, closed-shell UCJ-BCH2
optimizer described in :mod:`pyscf.ucj.ucj`.
"""

from pyscf.ucj.ucj import RUCJ, UCJ

__all__ = ['UCJ', 'RUCJ']
