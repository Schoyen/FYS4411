from libc.stdlib cimport malloc, free

cdef class Particles:
    cdef particles m_particles

    def __cinit__(self, int num_walkers):
        self.m_particles.particles = \
            <particle *> malloc(sizeof(particle)*num_walkers)

        if not self.m_particles.particles:
            raise Exception("particles-array was not allocated")
