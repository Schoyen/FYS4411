from coupled_cluster.schemes.ccd_sparse import CoupledClusterDoublesSparse
from coupled_cluster.schemes.ccd_optimized import CoupledClusterDoublesOptimized
from coupled_cluster.schemes.ccd import CoupledClusterDoubles
from coupled_cluster.hartree_fock.scf_rhf import scf_rhf
from coupled_cluster.matrix_elements.generate_matrices import (
    get_one_body_elements, get_coulomb_elements,
    get_antisymmetrized_elements, add_spin_to_one_body_elements,
    get_one_body_elements_spin
)
from coupled_cluster.matrix_elements.index_map import (
    generate_index_map, IndexMap
)
from coupled_cluster.hartree_fock.basis_transformation import (
    transform_one_body_elements, transform_two_body_elements
)

import numpy as np
import time
import os

ccd_method = CoupledClusterDoublesOptimized

file_path = os.path.join("..", "dat")
filename = os.path.join(file_path, "coulomb_{0}.pkl")

num_shells = 5
generate_index_map(num_shells)

omega = 1.0
l = IndexMap.shell_arr[-1]
n = 2
theta_hf = 0.1
theta_ho = 0.1
max_iterations = 1000
convergence = 1e-2

filename = filename.format(l)

print ("""
w = {0},
num_shells = {1},
l = {2},
n = {3},
theta_hf = {4},
theta_ho = {5},
filename = {6}
""".format(omega, num_shells, l, n, theta_hf, theta_ho, filename))

h = omega * get_one_body_elements(l)
t0 = time.time()
u = np.sqrt(omega) * get_coulomb_elements(l, filename=filename, tol=1e-12)
t1 = time.time()
print ("Time spent creating Coulomb elements: {0} sec".format(t1 - t0))

t0 = time.time()
c, energy = scf_rhf(h.todense(), u, np.eye(l//2), n//2, tol=convergence)
t1 = time.time()
print ("Time spent in SCF RHF: {0} sec".format(t1 - t0))

print ("\tRHF Energy: {0:.6f}".format(energy))

hi = transform_one_body_elements(h, c)
t0 = time.time()
oi = transform_two_body_elements(u, c)
t1 = time.time()
print ("Time spent transforming two body elements: {0} sec".format(t1 - t0))

_h = add_spin_to_one_body_elements(hi, l)
t0 = time.time()
_u = get_antisymmetrized_elements(l, oi=oi, tol=1e-12)
t1 = time.time()
print ("Time spent antisymmetrizing two body elements: {0} sec".format(t1 - t0))

t0 = time.time()
ccd_hf = ccd_method(_h.todense(), _u.todense(), n)
t1 = time.time()
print ("Time spent setting up {0} code with HF basis: {1} sec"
    .format(ccd_method.__name__, t1 - t0))

t0 = time.time()
energy, iterations = ccd_hf.compute_energy(
        tol=convergence, theta=theta_hf, max_iterations=max_iterations)
t1 = time.time()
print ("Time spent computing {0} energy with HF basis: {1:.4f} sec"
   .format(ccd_method.__name__, t1 - t0))
print (("\tCCD (HF) Energy: {0:.6f}\n\tIterations: {1}\n\tSecond/iteration:" \
        + "{2:.4f}").format(energy, iterations, (t1 - t0)/iterations))
