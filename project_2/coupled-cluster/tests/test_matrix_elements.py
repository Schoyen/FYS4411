import pytest
import sparse
from coupled_cluster.matrix_elements.index_map import (
        get_indices_nm, generate_index_map
)

from coupled_cluster.matrix_elements.coulomb_interface import (
        get_coulomb_element
)

from coupled_cluster.matrix_elements.generate_matrices import (
        get_coulomb_elements, get_antisymmetrized_elements,
        get_one_body_elements_spin
)

def test_two_body_generation_one():
    orbital_integrals = pytest.orbital_integrals
    l = pytest.l
    num_shells = pytest.num_shells

    generate_index_map(num_shells)

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

    #sparse.utils.assert_eq(
    #        orbital_integrals, get_coulomb_elements(l), atol=1e-5, rtol=1e-5)
    u = get_coulomb_elements(l)

    p, q, r, s = u.coords

    with open("tests/dat/foo.dat", "w") as f:
        for _p, _q, _r, _s in zip(p, q, r, s):
            f.write(
                    "{0} {1} {2} {3} {4}\n".format(
                        _p, _q, _r, _s, u[_p, _q, _r, _s]))

    #sparse.utils.assert_eq(orbital_integrals, u)

    for _p, _q, _r, _s in zip(p, q, r, s):
        assert abs(
                u[_p, _q, _r, _s] - orbital_integrals[_p, _q, _r, _s]) < 1e-5, \
                "(p, q, r, s) = ({0}, {1}, {2}, {3}), u = {4}, oi = {5}".format(
                        _p, _q, _r, _s, u[_p, _q, _r, _s],
                        orbital_integrals[_p, _q, _r, _s])
