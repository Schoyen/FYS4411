cdef extern from "wavefunction.h":
    cdef struct parameters:
        unsigned int num_parameters
        double parameters[]
