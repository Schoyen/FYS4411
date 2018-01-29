cdef extern from "wavefunction_config.h":
    cdef enum:
        DIMENSIONALITY

cdef extern from "particles.h":
    cdef struct particle:
        double position[DIMENSIONALITY]

    cdef struct particles:
        unsigned int num_particles
        particle *particles
