from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np

np.import_array()

cdef class Wavefunction:
    cdef wavefunction m_wavefunction

    def __cinit__(
            self, unsigned int num_particles, unsigned int dimensionality):

        self.m_wavefunction.num_particles = num_particles
        self.m_wavefunction.dimensionality = dimensionality

        allocate_variational_parameters(&self.m_wavefunction)
        allocate_particles(&self.m_wavefunction)

    def get_parameters(self):
        cdef unsigned int i, num_parameters

        num_parameters = self.m_wavefunction.parameters.num_parameters

        return [
            self.m_wavefunction.parameters.parameters[i]
            for i in range(num_parameters)]

    def set_parameters(self, list values):
        cdef unsigned int i, num_parameters

        num_parameters = self.m_wavefunction.parameters.num_parameters

        for i in range(num_parameters):
            self.m_wavefunction.parameters.parameters[i] = <double> values[i]


    def __dealloc__(self):
        free_variational_parameters(&self.m_wavefunction)
        free_particles(&self.m_wavefunction)

        return particles

    def initialize_particles(self, double spread):
        cdef unsigned int i

        for i in range(self.m_particles.num_particles):
            self.initialize_particle(&(self.m_particles.particles[i]), spread)

    cdef inline initialize_particle(self, particle *particle, double spread):
        cdef int i

        for i in range(DIMENSIONALITY):
            particle.position[i] = spread*(2.0*np.random.random() - 1.0)

    def __dealloc__(self):
        free(self.m_particles.particles)


def perform_metropolis(
        Wavefunction wavefunction, Particles particles, double step_length,
        unsigned int num_samples):

    cdef double energy

    energy = metropolis_sampling(
        &particles.m_particles,
        wavefunction.m_parameters,
        step_length,
        num_samples
    )

    return energy

def normalize_energies(
        np.ndarray[double, ndim=1] energies, unsigned int num_samples,
        unsigned int num_particles):

    return energies/(num_samples*num_particles)
