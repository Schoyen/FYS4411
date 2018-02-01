cdef extern from "wavefunction_config.h":
    cdef enum:
        DIMENSIONALITY

cdef extern from "metropolis_sampling.h":
    double perform_metropolis_step(
            particles *particles,
            parameters *parameters,
            double step_length)

    double metropolis_sampling(
            particles *particles,
            parameters *parameters,
            double step_length,
            unsigned int num_samples)

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

cdef extern from "particles.h":
    cdef struct particle:
        double position[DIMENSIONALITY]

    cdef struct particles:
        unsigned int num_particles
        particle *particles
