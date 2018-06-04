import numpy as np
import pytest

from coupled_cluster.schemes.ccd_sparse import CoupledClusterDoublesSparse
from coupled_cluster.schemes.ccd import CoupledClusterDoubles
from coupled_cluster.schemes.ccd_optimized import CoupledClusterDoublesOptimized
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

def test_compare_energy():
    l = pytest.l
    n = pytest.n
    num_shells = pytest.num_shells

    convergence_criteria = 1e-8

    generate_index_map(num_shells)
    h = get_one_body_elements(l)
    u = get_coulomb_elements(l)

    c, energy = scf_rhf(h.todense(), u, np.eye(l//2), n//2)

    hi = transform_one_body_elements(h, c)
    oi = transform_two_body_elements(u, c)

    _h = add_spin_to_one_body_elements(hi, l)
    _u = get_antisymmetrized_elements(l, oi=oi)

    ccd = CoupledClusterDoubles(_h.todense(), _u.todense(), n)
    ccd_sparse = CoupledClusterDoublesSparse(_h, _u, n)

    energy, iterations = ccd.compute_energy(tol=convergence_criteria)
    energy_sparse, iterations_sparse = ccd_sparse.compute_energy(
            tol=convergence_criteria)

    assert abs(energy_sparse - energy) < convergence_criteria

def test_optimized_energy():
    l = pytest.l
    n = pytest.n
    num_shells = pytest.num_shells

    convergence_criteria = 1e-8

    generate_index_map(num_shells)
    h = get_one_body_elements(l)
    u = get_coulomb_elements(l)

    c, energy = scf_rhf(h.todense(), u, np.eye(l//2), n//2)

    hi = transform_one_body_elements(h, c)
    oi = transform_two_body_elements(u, c)

    _h = add_spin_to_one_body_elements(hi, l).todense()
    _u = get_antisymmetrized_elements(l, oi=oi).todense()

    ccd = CoupledClusterDoubles(_h, _u, n)
    ccd_opt = CoupledClusterDoublesOptimized(_h, _u, n)

    energy, iterations = ccd.compute_energy(tol=convergence_criteria)
    energy_opt, iterations_opt = ccd_opt.compute_energy(
            tol=convergence_criteria)

    assert abs(energy - energy_opt) < convergence_criteria
