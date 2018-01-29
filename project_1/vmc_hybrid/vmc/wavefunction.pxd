from particles cimport DIMENSIONALITY

cdef extern from "wavefunction.h":
    cdef struct parameters:
        unsigned int num_parameters
        double *parameters

    parameters *get_variational_parameters()
    void free_parameters_struct(parameters *parameters)

    double local_energy(
            parameters *parameters, double position[DIMENSIONALITY])

    double ratio(
            parameters *parameters, double new_position[DIMENSIONALITY],
            double old_position[DIMENSIONALITY])

    double wavefunction(
            parameters *parameters, double position[DIMENSIONALITY])
