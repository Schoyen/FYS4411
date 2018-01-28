cdef extern from "particles.h":
    cdef enum:
        DIMENSIONALITY = 3

    cdef struct particle:
        double position[DIMENSIONALITY]

    cdef struct particles:
        particle *particles
