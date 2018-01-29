from particles cimport *
from wavefunction cimport *

cdef extern from "metropolis_sampling.h":
    double perform_metropolis_step(
            particles *particles,
            parameters *parameters,
            double step_length)

    double metropolis_sampling(
            particles *particles,
            parameters *parameters,
            double step_length,
            unsigned int num_samples)
