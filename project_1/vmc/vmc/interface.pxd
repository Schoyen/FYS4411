from libcpp cimport bool

cdef extern from "<valarray>" namespace "std":
    cdef cppclass valarray[T]:
        T& operator[](int)
        unsigned int size()

cdef extern from "constants.h":
    cdef double HBAR

cdef extern from "wavefunction.h":
    cdef cppclass Wavefunction:
        double get_mass()
        double get_frequency()
        unsigned int get_num_particles()
        unsigned int get_num_dimensions()
        double *get_parameters()

cdef extern from "simple_gaussian.h":
    cdef cppclass SimpleGaussian(Wavefunction):
        SimpleGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                double mass,
                double omega,
                double *parameters,
                double *particles) except +

cdef extern from "simple_gaussian_numerical.h":
    cdef cppclass SimpleGaussianNumerical(SimpleGaussian):
        SimpleGaussianNumerical(
                unsigned int num_particles,
                unsigned int num_dimensions,
                double mass,
                double omega,
                double h,
                double *parameters,
                double *particles) except +

cdef extern from "interacting_elliptical_gaussian.h":
    cdef cppclass InteractingEllipticalGaussian(Wavefunction):
        InteractingEllipticalGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                double mass,
                double omega,
                double beta,
                double radius,
                double *parameters,
                double *particles) except +

cdef extern from "hamiltonian.h":
    cdef cppclass Hamiltonian:
        pass

cdef extern from "harmonic_oscillator.h":
    cdef cppclass HarmonicOscillator(Hamiltonian):
        pass

cdef extern from "elliptical_harmonic_oscillator.h":
    cdef cppclass EllipticalHarmonicOscillator(Hamiltonian):
        EllipticalHarmonicOscillator(double _lambda) except +

cdef extern from "monte_carlo_method.h":
    cdef cppclass MonteCarloMethod:
        MonteCarloMethod() except +
        MonteCarloMethod(int seed) except +

        bool step(Wavefunction *wavefunction, double step_length)

cdef extern from "sampler.h":
    cdef cppclass Sampler:
        Sampler(
                Wavefunction *wavefunction,
                Hamiltonian *hamiltonian,
                MonteCarloMethod *solver) except +

        void sample(unsigned int num_samples, double step_length,
                double *local_energies)
        void set_one_body_parameters(
                double r_min, double r_max, unsigned int num_bins, double *bins)
        double get_variance()
        double get_energy()
        double get_energy_squared()
        double get_acceptance_ratio()
        valarray[double] get_variational_parameters_gradient()

cdef extern from "metropolis_algorithm.h":
    cdef cppclass MetropolisAlgorithm(MonteCarloMethod):
        MetropolisAlgorithm() except +
        MetropolisAlgorithm(int seed) except +

cdef extern from "importance_metropolis.h":
    cdef cppclass ImportanceMetropolis(MonteCarloMethod):
        ImportanceMetropolis(double diffusion_coefficient) except +
        ImportanceMetropolis(double diffusion_coefficient, int seed) except +
