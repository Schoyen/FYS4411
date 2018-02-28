import numpy as np
cimport numpy as np

np.import_array()

cdef class PyWavefunction:

    cdef Wavefunction *wavefunction
    cdef double[:, ::1] particles
    cdef double[::1] parameters
    cdef unsigned int num_particles
    cdef unsigned int num_dimensions
    cdef unsigned int num_parameters
    cdef double spread

    def __init__(self, unsigned int num_particles,
            unsigned int num_dimensions, unsigned int num_parameters,
            double spread=1.0):

        self.num_particles = num_particles
        self.num_dimensions = num_dimensions
        self.num_parameters = num_parameters
        self.spread = spread

        self.particles = spread \
                * (2*np.random.random((num_particles, num_dimensions)) - 1.0)
        self.parameters = np.zeros(num_parameters)

    def redistribute(self, double spread=-1):

        cdef np.ndarray[double, ndim=2, mode="c"] distro
        cdef unsigned int i, j

        if spread < 0:
            spread = self.spread

        distro = spread \
                * (2*np.random.random(
                    (self.num_particles, self.num_dimensions)) - 1.0)

        for i in range(self.num_particles):
            for j in range(self.num_dimensions):
                self.particles[i, j] = distro[i, j]

    def get_particles(self):
        return np.asarray(self.particles)

    def set_parameters(self, np.ndarray[double, ndim=1, mode="c"] parameters):
        cdef unsigned int i

        if parameters.size != self.num_parameters:
            raise Exception(
                    "The length of the parameters array must be equal to "
                    + "the number of parameters")

        for i in range(self.num_parameters):
            self.parameters[i] = parameters[i]

    def get_parameters(self):
        return np.asarray(self.parameters)

    def evaluate(self):
        return self.wavefunction.evaluate()

    def compute_laplacian(self):
        return self.wavefunction.compute_laplacian()

    def compute_position_squared_sum(self):
        return self.wavefunction.compute_position_squared_sum()

cdef class PySimpleGaussian(PyWavefunction):

    def __init__(self, unsigned int num_particles, unsigned int num_dimensions,
            unsigned int num_parameters, double mass, double omega,
            double spread=1.0):

        super().__init__(num_particles, num_dimensions, num_parameters, spread)

        self.wavefunction = new SimpleGaussian(
                num_particles, num_dimensions, num_parameters, mass, omega,
                &self.parameters[0], &self.particles[0, 0])

cdef class PySimpleGaussianNumerical(PyWavefunction):

    def __init__(self, unsigned int num_particles, unsigned int num_dimensions,
            unsigned int num_parameters, double mass, double omega,
            double h=1e-7, double spread=1.0):
        super().__init__(num_particles, num_dimensions, num_parameters, spread)

        self.wavefunction = new SimpleGaussianNumerical(
                num_particles, num_dimensions, num_parameters, mass, omega, h,
                &self.parameters[0], &self.particles[0, 0])

cdef class PyHamiltonian:
    cdef Hamiltonian *hamiltonian

    def compute_local_energy(self, PyWavefunction wavefunction):
        return self.hamiltonian.compute_local_energy(
                wavefunction.wavefunction)

    def compute_potential_energy(self, PyWavefunction wavefunction):
        return self.hamiltonian.compute_potential_energy(
                wavefunction.wavefunction)

cdef class PyHarmonicOscillator(PyHamiltonian):

    def __cinit__(self, double mass, double omega):
        self.hamiltonian = new HarmonicOscillator(mass, omega)

cdef class PyMonteCarloMethod:
    cdef MonteCarloMethod *method

cdef class PySampler:
    cdef Sampler *sampler
    cdef double[::1] local_energies

    def __cinit__(self, PyWavefunction wavefunction, PyHamiltonian hamiltonian,
            PyMonteCarloMethod solver, unsigned int num_local_energies,
            unsigned int stride_local_energies):

        self.wavefunction = wavefunction
        self.hamiltonian = hamiltonian
        self.solver = solver
        self.local_energies = np.zeros(num_local_energies)

        self.sampler = new Sampler(
                wavefunction.wavefunction,
                hamiltonian.hamiltonian,
                solver.method,
                num_local_energies,
                stride_local_energies,
                &self.local_energies[0])

    def get_local_energies(self):
        return np.asarray(self.local_energies)

    def get_variance(self):
        return self.sampler.get_variance()

    def get_energy(self):
        return self.sampler.get_energy()

    def get_energy_squared(self):
        return self.sampler.get_energy_squared()

    def get_ratio_of_accepted_steps(self):
        return self.sampler.get_ratio_of_accepted_steps()

cdef class PyMetropolisAlgorithm(PyMonteCarloMethod):

    def __cinit__(self, unsigned int num_particles):
        self.method = new MetropolisAlgorithm(num_particles)

    #def step(self, PyWavefunction wavefunction, double step_length):
    #    return self.method.step(wavefunction.wavefunction, step_length)

    #def run(self, PyWavefunction wavefunction, PyHamiltonian hamiltonian,
    #        double step_length, unsigned int num_samples):
    #    return self.method.run(
    #            wavefunction.wavefunction, hamiltonian.hamiltonian,
    #            step_length, num_samples)

    #def run_with_variance(self, PyWavefunction wavefunction,
    #        PyHamiltonian hamiltonian, double step_length,
    #        unsigned int num_samples):
    #    cdef double variance
    #    cdef double energy

    #    variance = 0
    #    energy = self.method.run_variance(
    #            wavefunction.wavefunction, hamiltonian.hamiltonian,
    #            step_length, num_samples, &variance)

    #    return energy, variance

    #def sample_local_energy(self, PyWavefunction wavefunction,
    #        PyHamiltonian hamiltonian, double step_length,
    #        unsigned int num_samples):

    #    cdef double[::1] local_energies
    #    cdef double energy

    #    local_energies = np.zeros(num_samples)

    #    energy = self.method.run(
    #            wavefunction.wavefunction, hamiltonian.hamiltonian,
    #            step_length, num_samples, &local_energies[0])

    #    return energy, np.asarray(local_energies)

cdef class PySteepestDescentMetropolis(PyMonteCarloMethod):

    def __cinit__(self, unsigned int num_particles):
        self.method = new SteepestDescentMetropolis(num_particles)

    #def steepest_descent(self, PyWavefunction wavefunction,
    #        PyHamiltonian hamiltonian, double gamma, unsigned int num_samples):
    #    return self.method.steepest_descent(
    #            wavefunction.wavefunction, hamiltonian.hamiltonian, gamma,
    #            num_samples)

cdef class PyImportanceMetropolis(PyMonteCarloMethod):

    def __cinit__(self, unsigned int num_particles):
        self.method = new ImportanceMetropolis(num_particles)

#    def step(self, PyWavefunction wavefunction, double step_length):
#        return self.method.step(wavefunction.wavefunction, step_length)
#
    #def run(self, PyWavefunction wavefunction, PyHamiltonian hamiltonian,
    #        double step_length, unsigned int num_samples):
    #    return self.method.run(
    #            wavefunction.wavefunction, hamiltonian.hamiltonian,
    #            step_length, num_samples)
