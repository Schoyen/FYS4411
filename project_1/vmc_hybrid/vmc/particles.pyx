from libc.stdlib cimport malloc, free

cdef class Particles:
    cdef particles *m_particles

    def __cinit__(self, int num_walkers):
        m_particles = <particles *> malloc(sizeof(particles))
        m_particles.particles = \
            <particle *> malloc(sizeof(particle)*num_walkers)
