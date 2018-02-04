cdef extern from "wavefunction.h":
    cdef struct particle:
        double *position

    cdef struct particles:
        particle_t *particles

    cdef struct parameters:
        unsigned int num_parameters
        double *parameters

    cdef struct wavefunction:
        unsigned int num_particles
        unsigned int dimensionality
        parameters *parameters
        particles *particles

    void allocate_variational_parameters(wavefunction *wavefunction)
    void free_variational_parameters(wavefunction *wavefunction)

    void allocate_particles(wavefunction *wavefunction)
    void free_particles(wavefunction *wavefunction)

cdef extern from "metropolis_sampling.h":
    double perform_metropolis_step(
            wavefunction *wavefunction, double step_length)

    double metropolis_sampling(
            wavefunction *wavefunction, double step_length,
            unsigned int num_samples)
