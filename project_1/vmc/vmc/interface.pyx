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

    def set_parameters(self, np.ndarray[double, ndim=1, mode="c"] parameters):
        cdef unsigned int i

        if parameters.size != self.num_parameters:
            raise Exception(
                    "The length of the parameters array must be equal to "
                    + "the number of parameters")

        for i in range(self.num_parameters):
            self.parameters[i] = parameters[i]

    def get_mass(self):
        return self.wavefunction.get_mass()

    def get_frequency(self):
        return self.wavefunction.get_frequency()

    def get_num_particles(self):
        return self.wavefunction.get_num_particles()

    def get_num_dimensions(self):
        return self.wavefunction.get_num_dimensions()

    def get_parameters(self):
        return np.asarray(self.parameters)

    def __dealloc__(self):
        del self.wavefunction

cdef class PySimpleGaussian(PyWavefunction):

    def __init__(self, unsigned int num_particles, unsigned int num_dimensions,
            double mass, double omega, double spread=1.0):

        cdef unsigned int num_parameters = 1
        super().__init__(num_particles, num_dimensions, num_parameters, spread)

        self.wavefunction = new SimpleGaussian(
                num_particles, num_dimensions, mass, omega,
                &self.parameters[0], &self.particles[0, 0])

cdef class PySimpleGaussianNumerical(PyWavefunction):

    def __init__(self, unsigned int num_particles, unsigned int num_dimensions,
            double mass, double omega, double h=1e-7, double spread=1.0):

        cdef unsigned int num_parameters = 1
        super().__init__(num_particles, num_dimensions, num_parameters, spread)

        self.wavefunction = new SimpleGaussianNumerical(
                num_particles, num_dimensions, mass, omega, h,
                &self.parameters[0], &self.particles[0, 0])

cdef class PyInteractingEllipticalGaussian(PyWavefunction):
    cdef double radius

    def __init__(self, unsigned int num_particles, unsigned int num_dimensions,
            double mass, double omega, double beta, double radius,
            double spread=1.0):

        assert radius < spread, "Radius must be smaller than spread"

        cdef unsigned int num_parameters = 1

        self.radius = radius
        super().__init__(num_particles, num_dimensions, num_parameters, spread)

        self._remove_overlap()

        self.wavefunction = new InteractingEllipticalGaussian(
                num_particles, num_dimensions, mass, omega, beta, radius,
                &self.parameters[0], &self.particles[0, 0])

    def _remove_overlap(self):
        cdef unsigned int p_i, i, p_j, j

        for p_i in range(self.num_particles):
            for p_j in range(self.num_particles):
                if p_i == p_j:
                    continue

                for i in range(self.num_dimensions):

                    while abs(self.particles[p_i, i] - self.particles[p_j, i]) \
                            <= self.radius:
                        self.particles[p_j, i] = \
                                self.spread*(2.0*np.random.random() - 1.0)

    def redistribute(self, double spread=-1):
        super().redistribute(spread=spread)
        self._remove_overlap()

cdef class PyHamiltonian:
    cdef Hamiltonian *hamiltonian

cdef class PyHarmonicOscillator(PyHamiltonian):

    def __cinit__(self):
        self.hamiltonian = new HarmonicOscillator()

    def compute_exact_energy(self, PyWavefunction wavefunction, alphas=[]):
        cdef double mass, omega
        cdef unsigned int num_dimensions, num_particles
        cdef np.ndarray[double, ndim=1, mode="c"] alpha

        if len(alphas) == 0:
            alpha = wavefunction.get_parameters()
        else:
            alpha = alphas

        mass = wavefunction.get_mass()
        omega = wavefunction.get_frequency()
        num_dimensions = wavefunction.get_num_dimensions()
        num_particles = wavefunction.get_num_particles()

        energy = HBAR**2*alpha/(2.0*mass) + mass*omega**2/(8.0*alpha)
        energy *= num_dimensions*num_particles

        return energy

cdef class PyEllipticalHarmonicOscillator(PyHamiltonian):

    def __cinit__(self, double _lambda):
        self.hamiltonian = new EllipticalHarmonicOscillator(_lambda)

cdef class PyMonteCarloMethod:
    cdef MonteCarloMethod *method

cdef class PySampler:
    cdef Sampler *sampler
    cdef double[::1] local_energies
    cdef double[::1] bins
    cdef PyWavefunction wavefunction

    def __init__(self, PyWavefunction wavefunction, PyHamiltonian hamiltonian,
            PyMonteCarloMethod solver):

        self.local_energies = np.zeros(1)
        self.wavefunction = wavefunction

        self.sampler = new Sampler(
                wavefunction.wavefunction,
                hamiltonian.hamiltonian,
                solver.method)

    def initialize(self, np.ndarray[double, ndim=1, mode="c"] parameters):
        self.wavefunction.redistribute()
        self.wavefunction.set_parameters(parameters)

    def get_local_energies(self):
        return np.asarray(self.local_energies)

    def get_variance(self):
        return self.sampler.get_variance()

    def get_stddev(self):
        return np.sqrt(self.get_variance())

    def get_energy(self):
        return self.sampler.get_energy()

    def get_energy_squared(self):
        return self.sampler.get_energy_squared()

    def get_acceptance_ratio(self):
        return self.sampler.get_acceptance_ratio()

    def get_parameter_gradient(self):
        cdef np.ndarray[double, ndim=1, mode="c"] gradient
        cdef valarray[double] _gradient
        cdef unsigned int i

        _gradient = self.sampler.get_variational_parameters_gradient()
        gradient = np.zeros(_gradient.size())

        for i in range(gradient.size):
            gradient[i] = _gradient[i]

        return gradient

    def get_one_body_densities(self):
        return np.asarray(self.bins)

    def initialize_one_body_densities(self, double r_min, double r_max,
            unsigned int num_bins):
        self.bins = np.zeros(num_bins)

        self.sampler.set_one_body_parameters(
                r_min, r_max, num_bins, &self.bins[0])

    def find_minimum(self,
            np.ndarray[double, ndim=1, mode="c"] start_parameters,
            double gamma, double gamma_scale, unsigned int max_iterations,
            *sample_args, quiet=False, double magnitude_limit=5,
            double tol=1e-9, **sample_kwargs):

        cdef np.ndarray[double, ndim=1, mode="c"] parameters, \
                parameters_prev, gradient, gradient_prev
        cdef unsigned int i

        parameters = start_parameters
        breadcrumbs = [parameters.copy()]

        gradient_prev = np.ones(parameters.size)*1e10

        for i in range(max_iterations):
            self.initialize(parameters)
            self.sample(*sample_args, **sample_kwargs)

            gradient = self.get_parameter_gradient()

            if not quiet:
                print(
"""
Iteration: {0}
Gradient: {1}
Old parameters: {2}""".format(i, gradient, parameters))

            if np.all(gradient**2 > gradient_prev**2):
                gamma *= gamma_scale

            parameters_prev = parameters.copy()
            parameters += -gamma*gradient

            mask = np.abs(parameters) > magnitude_limit*np.abs(parameters_prev)
            while np.any(mask):
                parameters[mask] /= magnitude_limit
                mask = np.abs(parameters) \
                        > magnitude_limit*np.abs(parameters_prev)

            # We only use positive variational parameters
            mask = parameters < 0
            parameters[mask] = -parameters[mask]
            breadcrumbs.append(parameters.copy())

            if not quiet:
                print("""Current parameters: {0}""".format(parameters))

            if np.all((parameters - parameters_prev)**2 < tol):
                break

        return np.asarray(breadcrumbs)

    def sample(self, unsigned int num_samples, double step_length,
            unsigned int num_thermalization_steps=0,
            sample_local_energies=False):

        if num_thermalization_steps > 0:
            self.sampler.sample(
                    num_thermalization_steps, step_length, <double *> 0)

        if sample_local_energies:
            if len(self.local_energies) != num_samples:
                self.local_energies = np.zeros(num_samples)

            self.sampler.sample(
                    num_samples, step_length, &self.local_energies[0])
        else:
            self.sampler.sample(num_samples, step_length, <double *> 0)

cdef class PyMetropolisAlgorithm(PyMonteCarloMethod):

    def __cinit__(self, seed=None):
        if seed:
            self.method = new MetropolisAlgorithm(seed)
        else:
            self.method = new MetropolisAlgorithm()

cdef class PyImportanceMetropolis(PyMonteCarloMethod):

    def __cinit__(self, double diffusion_coefficient=0.5, seed=None):
        if seed:
            self.method = new ImportanceMetropolis(diffusion_coefficient, seed)
        else:
            self.method = new ImportanceMetropolis(diffusion_coefficient)
