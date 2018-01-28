from libc.stdlib cimport malloc, free

cdef class Particles:
    cdef particles m_particles

    def __cinit__(self, unsigned int num_particles):
        self.m_particles.num_particles = num_particles
        self.m_particles.particles = \
            <particle *> malloc(sizeof(particle)*num_particles)

        if not self.m_particles.particles:
            raise Exception("particles-array was not allocated")

    def __dealloc__(self):
        free(self.m_particles.particles)
