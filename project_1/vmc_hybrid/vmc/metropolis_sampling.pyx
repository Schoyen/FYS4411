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

    def initialize_wavefunction(self, double spread):
        # TODO: Distribute particles
        # TODO: Evaluate wavefunction and store in last_value
        pass

    def get_parameters(self):
        cdef unsigned int i, num_parameters

        num_parameters = self.m_wavefunction.num_parameters

        return [
            self.m_wavefunction.parameters[i]
            for i in range(num_parameters)]

    def set_parameters(self, list values):
        cdef unsigned int i, num_parameters

        num_parameters = self.m_wavefunction.num_parameters

        for i in range(num_parameters):
            self.m_wavefunction.parameters[i] = <double> values[i]


    def __dealloc__(self):
        free_variational_parameters(&self.m_wavefunction)
        free_particles(&self.m_wavefunction)

def perform_varying_metropolis(
        Wavefunction wavefunction, double step_length,
        unsigned int num_samples, np.ndarray[double, ndim=2] parameters):

    cdef np.ndarray[double, ndim=1] energies
    cdef unsigned int i

    energies = np.zeros(len(parameters))

    for i in range(len(parameters)):
        wavefunction.set_parameters(parameters[i])
        energies[i] = _perform_metropolis(
                wavefunction, step_length, num_samples)

    return energies

cdef double _perform_metropolis(
        Wavefunction wavefunction, double step_length,
        unsigned int num_samples):

    return metropolis_sampling(
            &wavefunction.m_wavefunction, step_length, num_samples)

def perform_metropolis(
        Wavefunction wavefunction, double step_length,
        unsigned int num_samples):

    return _perform_metropolis(wavefunction, step_length, num_samples)

def normalize_energies(
        np.ndarray[double, ndim=1] energies, unsigned int num_samples,
        unsigned int num_particles):

    return energies/(num_samples*num_particles)
