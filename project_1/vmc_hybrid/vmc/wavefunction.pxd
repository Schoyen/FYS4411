cdef extern from "wavefunction.h":
    cdef struct parameters:
        unsigned int num_parameters
        double *parameters

    parameters *get_variational_parameters()
    void free_parameters_struct(parameters *parameters)
