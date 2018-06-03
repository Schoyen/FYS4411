import numpy as np
import pytest

from coupled_cluster.schemes.ccd_sparse import CoupledClusterDoublesSparse
from coupled_cluster.schemes.ccd import CoupledClusterDoubles
from coupled_cluster.hartree_fock.scf_rhf import scf_rhf
from coupled_cluster.matrix_elements.generate_matrices import (
        get_one_body_elements, get_coulomb_elements,
        get_antisymmetrized_elements, add_spin_to_one_body_elements
)
from coupled_cluster.hartree_fock.basis_transformation import (
        transform_one_body_elements, transform_two_body_elements
)
from coupled_cluster.matrix_elements.index_map import (
        generate_index_map
)

def test_hartree_fock_energy_sparse():
    l = pytest.l
    n = pytest.n
    num_shells = pytest.num_shells

    generate_index_map(num_shells)
    h = get_one_body_elements(l)
    u = get_coulomb_elements(l)

    c, energy = scf_rhf(h.todense(), u, np.eye(l//2), n//2)

    hi = transform_one_body_elements(h, c)
    oi = transform_two_body_elements(u, c)

    _h = add_spin_to_one_body_elements(hi, l)
    _u = get_antisymmetrized_elements(l, oi=oi)

    ccd = CoupledClusterDoublesSparse(_h, _u, n)
    assert abs(ccd.compute_reference_energy() - energy) < 1e-8

def test_hartree_fock_energy():
    l = pytest.l
    n = pytest.n
    num_shells = pytest.num_shells

    generate_index_map(num_shells)
    h = get_one_body_elements(l)
    u = get_coulomb_elements(l)

    c, energy = scf_rhf(h.todense(), u, np.eye(l//2), n//2)

    hi = transform_one_body_elements(h, c)
    oi = transform_two_body_elements(u, c)

    _h = add_spin_to_one_body_elements(hi, l)
    _u = get_antisymmetrized_elements(l, oi=oi)

    ccd = CoupledClusterDoubles(_h.todense(), _u.todense(), n)
    assert abs(ccd.compute_reference_energy() - energy) < 1e-8
