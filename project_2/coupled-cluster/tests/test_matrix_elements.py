import pytest
import sparse
from coupled_cluster.matrix_elements.index_map import (
        get_indices_nm
)

from coupled_cluster.matrix_elements.coulomb_interface import (
        get_coulomb_element
)

from coupled_cluster.matrix_elements.generate_matrices import (
        get_coulomb_elements, get_antisymmetrized_elements,
        get_one_body_elements
)

def test_two_body_generation_one():
    orbital_integrals = pytest.orbital_integrals
    l = pytest.l

    _p, _q, _r, _s = orbital_integrals.coords

    for p in range(l//2):
        for q in range(l//2):
            for r in range(l//2):
                for s in range(l//2):
                    gen_val = get_coulomb_element(
                            *get_indices_nm(p),
                            *get_indices_nm(q),
                            *get_indices_nm(r),
                            *get_indices_nm(s)
                    )

                    assert abs(orbital_integrals[p, q, r, s] - gen_val) < 1e-5

def test_two_body_generation():
    orbital_integrals = pytest.orbital_integrals
    l = pytest.l

    sparse.utils.assert_eq(
            orbital_integrals, get_coulomb_elements(l), atol=1e-5, rtol=1e-5)

def test_two_body_antisymmetric_generation():
    u = pytest.u
    l = pytest.l

    sparse.utils.assert_eq(
            u, get_antisymmetrized_elements(l), atol=1e-5, rtol=1e-5)

def test_one_body_generation():
    h = pytest.h
    l = pytest.l

    sparse.utils.assert_eq(
            h, get_one_body_elements(l), atol=1e-5, rtol=1e-5)
