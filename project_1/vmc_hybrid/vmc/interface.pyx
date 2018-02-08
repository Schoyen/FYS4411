import numpy as np
cimport numpy as np

np.import_array()

cdef class PyWavefunction:

    cdef Wavefunction *wavefunction


    def __cinit__(self, unsigned int num_particles,
            unsigned int num_dimensions, unsigned int num_parameters):
        cdef np.ndarray[double, ndim=2, mode="c"] self.particles
        cdef np.ndarray[double, ndim=1, mode="c"] self.parameters

        self.particles = np.random.random((num_particles, num_dimensions))
        self.parameters = np.zeros(num_parameters)

cdef class PySimpleGaussian(PyWavefunction):

    def __cinit__(self, unsigned int num_particles,
            unsigned int num_dimensions, unsigned int num_parameters):

        super.__cinit__(num_particles, num_dimensions, num_parameters)

        self.wavefunction = new SimpleGaussian(
                num_particles, num_dimensions, num_parameters,
                &self.parameters[0], &self.parameters[0, 0])

    def __dealloc__(self):
        del self.wavefunction
