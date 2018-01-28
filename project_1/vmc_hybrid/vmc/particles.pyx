from libc.stdlib cimport malloc, free

cdef class Particles:
    cdef particles *m_particles

    def __cinit__(self):
        m_particles = NULL
