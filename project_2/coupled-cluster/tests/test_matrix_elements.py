import pytest
from coupled_cluster.matrix_elements.coulomb_interface import (
    get_coulomb_element
)

def test_two_body_generation_one():
    orbital_integrals = pytest.orbital_integrals

    _p, _q, _r, _s = orbital_integrals.coords

    for p, q, r, s in zip(_p, _q, _r, _s):
        pass
