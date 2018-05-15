import numpy as np
import pytest

from coupled_cluster.schemes.ccd_sparse import CoupledClusterDoublesSparse
from coupled_cluster.hartree_fock.scf_rhf import scf_rhf
from coupled_cluster.matrix_elements.basis_transformation import (
        transform_one_body_elements, transform_two_body_elements
)

def test_hartree_fock_energy():
    h = pytest.h
    u = pytest.u
    l = pytest.l
    n = pytest.n

    c, energy = scf_rhf(h.todense(), u, np.eye(l), l, n)
    _h = transform_one_body_elements(h, c)
    print (_h.todense())

    _u = transform_two_body_elements(u, c)

    print (energy)

    ccd = CoupledClusterDoublesSparse(h, u, n)

    print (ccd.compute_reference_energy())

    ccd = CoupledClusterDoublesSparse(_h, _u, n)

    print (ccd.compute_reference_energy())

    wat
