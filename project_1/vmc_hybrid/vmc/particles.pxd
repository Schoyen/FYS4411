cdef extern from "particles.h":
    cdef struct particle:
        double *position

    cdef struct particles:
        particle *particles
