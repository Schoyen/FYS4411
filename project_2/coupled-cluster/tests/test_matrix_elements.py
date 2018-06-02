import pytest
import sparse
from coupled_cluster.matrix_elements.index_map import (
        get_indices_nm, generate_index_map
)

from coupled_cluster.matrix_elements.generate_matrices import (
        get_coulomb_elements, get_antisymmetrized_elements,
        get_one_body_elements_spin
)

def test_two_body_generation():
    orbital_integrals = pytest.orbital_integrals
    l = pytest.l
    num_shells = pytest.num_shells

    generate_index_map(num_shells)

    sparse.utils.assert_eq(
            orbital_integrals, get_coulomb_elements(l), atol=1e-5, rtol=1e-5)

def test_two_body_antisymmetric_generation():
    u = pytest.u
    l = pytest.l
    num_shells = pytest.num_shells

    generate_index_map(num_shells)

    sparse.utils.assert_eq(
            u, get_antisymmetrized_elements(l), atol=1e-5, rtol=1e-5)

def test_one_body_generation():
    h = pytest.h
    l = pytest.l
    num_shells = pytest.num_shells

    generate_index_map(num_shells)

    sparse.utils.assert_eq(
            h, get_one_body_elements_spin(l), atol=1e-5, rtol=1e-5)

def test_large_file():
    l = pytest.large_l
    orbital_integrals = pytest.large_oi
    num_shells = pytest.large_num_shells

    generate_index_map(num_shells)

    u = get_coulomb_elements(l)
    sparse.utils.assert_eq(orbital_integrals, u)
