import pytest
from coupled_cluster.matrix_elements.index_map import (
        get_indices_nm
)

from coupled_cluster.matrix_elements.coulomb_interface import (
        get_coulomb_element
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
