import numpy as np
cimport numpy as np

np.import_array()

cdef class PyWavefunction:

    cdef Wavefunction *wavefunction
    cdef double[:, ::1] particles
    cdef double[:] parameters
    cdef unsigned int num_particles
    cdef unsigned int num_dimensions
    cdef unsigned int num_parameters


    def __cinit__(self, unsigned int num_particles,
            unsigned int num_dimensions, unsigned int num_parameters):

        self.num_particles = num_particles
        self.num_dimensions = num_dimensions
        self.num_parameters = num_parameters

        self.particles = np.random.random((num_particles, num_dimensions))
        self.parameters = np.zeros(num_parameters)

    def get_particles(self):
        return np.asarray(self.particles)

    def get_parameters(self):
        return np.asarray(self.parameters)


cdef class PySimpleGaussian(PyWavefunction):

    def __cinit__(self, unsigned int num_particles,
            unsigned int num_dimensions, unsigned int num_parameters):

        self.wavefunction = new SimpleGaussian(
                num_particles, num_dimensions, num_parameters,
                &self.parameters[0], &self.particles[0, 0])
