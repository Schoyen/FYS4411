from coupled_cluster.schemes.ccd_sparse import CoupledClusterDoublesSparse
from coupled_cluster.schemes.ccd_optimized import CoupledClusterDoublesOptimized
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

file_path = os.path.join("..", "dat")
filename = os.path.join(file_path, "coulomb_{0}.pkl")

num_shells = 12
generate_index_map(num_shells)

omega = 1.0
l = IndexMap.shell_arr[-1]
n = 6
theta = 0.3

filename = filename.format(l)

print ("""
w = {0},
num_shells = {1},
l = {2},
n = {3},
theta = {4},
filename = {5}
""".format(omega, num_shells, l, n, theta, filename))

h = omega * get_one_body_elements(l)
t0 = time.time()
u = np.sqrt(omega) * get_coulomb_elements(l, filename=filename, tol=1e-12)
t1 = time.time()
print ("Time spent creating Coulomb elements: {0} sec".format(t1 - t0))

t0 = time.time()
c, energy = scf_rhf(h.todense(), u, np.eye(l//2), n//2)
t1 = time.time()
print ("Time spent in SCF RHF: {0} sec".format(t1 - t0))

print ("\tRHF Energy: {0}".format(energy))

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
ccd_hf_sparse = CoupledClusterDoublesSparse(_h, _u, n)
t1 = time.time()
print ("Time spent setting up CCD code with HF basis: {0} sec".format(t1 - t0))

t0 = time.time()
energy, iterations = ccd_hf_sparse.compute_energy(tol=1e-4, theta=theta)
t1 = time.time()
print ("Time spent computing CCD energy with HF basis: {0} sec".format(t1 - t0))
print ("\tCCD (HF) Energy: {0}\n\tIterations: {1}\n\tSecond/iteration: {2}"
    .format(energy, iterations, (t1 - t0)/iterations))

t0 = time.time()
ccd_hf = CoupledClusterDoublesOptimized(
        _h.todense(), _u.todense(), n, parallel=True)
t1 = time.time()
print ("Time spent setting up CCD (opt, parallel) code with HF basis: {0} sec"
    .format(t1 - t0))

t0 = time.time()
energy, iterations = ccd_hf.compute_energy(tol=1e-4, theta=theta)
t1 = time.time()
print ("Time spent computing CCD (opt, parallel) energy with HF basis: {0} sec"
   .format(t1 - t0))
print ("\tCCD (HF) Energy: {0}\n\tIterations: {1}\n\tSecond/iteration: {2}"
    .format(energy, iterations, (t1 - t0)/iterations))

__h = omega * get_one_body_elements_spin(l)
t0 = time.time()
__u = np.sqrt(omega) * get_antisymmetrized_elements(l, filename=filename)
t1 = time.time()
print ("Time spent getting antisymmetric two body elements: {0} sec".format(t1 - t0))

t0 = time.time()
ccd = CoupledClusterDoublesSparse(__h, __u, n)
t1 = time.time()
print ("Time spent setting up CCD code with HO basis: {0} sec".format(t1 - t0))
t0 = time.time()
energy, iterations = ccd.compute_energy(theta=theta)
t1 = time.time()
print ("Time spent computing CCD energy with HO basis: {0} sec".format(t1 - t0))
print ("\tCCD Energy: {0}\n\tIterations: {1}\n\tSecond/iteration: {2}"
    .format(energy, iterations, (t1 - t0)/iterations))

print (h.density)
print (u.density)
print (_u.density)
print (_h.density)
print (__u.density)
print (__h.density)
