from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np

cdef class Particles:
    cdef particles m_particles

    def __cinit__(self, unsigned int num_particles):
        self.m_particles.num_particles = num_particles
        self.m_particles.particles = \
            <particle *> malloc(sizeof(particle)*num_particles)

        if not self.m_particles.particles:
            raise Exception("particles-array was not allocated")

    def initialize_particles(self, double spread):
        cdef unsigned int i

        for i in range(self.m_particles.num_particles):
            self.initialize_particles(self.m_particles.particles[i], spread)

    cdef inline intialize_particle(self, particle *particle, double spread):
        cdef int i

        for i in range(DIMENSIONALITY):
            particle.position[i] = spread*(2.0*np.random.random() - 1.0)

    def __dealloc__(self):
        free(self.m_particles.particles)
