cdef extern from "wavefunction.h":
    cdef struct wavefunction:
        unsigned int num_particles
        unsigned int dimensionality
        unsigned int num_parameters

        double last_value

        double *parameters
        double **particles

    void allocate_variational_parameters(wavefunction *wavefunction)
    void free_variational_parameters(wavefunction *wavefunction)

    void allocate_particles(wavefunction *wavefunction)
    void free_particles(wavefunction *wavefunction)

    double evaluate_wavefunction(wavefunction *wavefunction)

cdef extern from "metropolis_sampling.h":
    double perform_metropolis_step(
            wavefunction *wavefunction, double step_length)

    double metropolis_sampling(
            wavefunction *wavefunction, double step_length,
            unsigned int num_samples)
