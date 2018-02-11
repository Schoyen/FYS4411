from libcpp cimport bool

cdef extern from "constants.h":
    cdef unsigned int HBAR

cdef extern from "wavefunction.h":
    cdef cppclass Wavefunction:
        Wavefunction(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double *parameters,
                double *particles) except +

        double compute_position_squared_sum()
        unsigned int get_num_dimensions()

        double evaluate()
        double compute_laplacian()

        void reset_particle_position(
                double position[], unsigned int particle_index)

        void copy_particle_position(
                double position[], unsigned int particle_index)

        void add_step(
                double step, unsigned int particle_index,
                unsigned int coordinate)

cdef extern from "simple_gaussian.h":
    cdef cppclass SimpleGaussian(Wavefunction):
        SimpleGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double *parameters,
                double *particles) except +

cdef extern from "hamiltonian.h":
    cdef cppclass Hamiltonian:
        pass

cdef extern from "harmonic_oscillator.h":
    cdef cppclass HarmonicOscillator(Hamiltonian):
        HarmonicOscillator(double mass, double omega) except +
        double compute_local_energy(Wavefunction *wavefunction)
        double compute_potential_energy(Wavefunction *wavefunction)

cdef extern from "monte_carlo_method.h":
    cdef cppclass MonteCarloMethod:
        MonteCarloMethod(unsigned int num_particles) except +

cdef extern from "metropolis_algorithm.h":
    cdef cppclass MetropolisAlgorithm(MonteCarloMethod):
        bool step(Wavefunction *wavefunction, double step_length)

        double run(
                Wavefunction *wavefunction, Hamiltonian *hamiltonian,
                double step_length, unsigned int num_samples)
