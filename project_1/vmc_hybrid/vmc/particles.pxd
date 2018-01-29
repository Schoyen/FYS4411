cdef extern from "particles.h":
    cdef enum:
        DIMENSIONALITY

    cdef struct particle:
        double position[DIMENSIONALITY]

    cdef struct particles:
        unsigned int num_particles
        particle *particles
