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

    def __cinit__(self, unsigned int num_particles,
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
        
        if spread < 0:
            spread = self.spread

        self.particles = spread \
             * (2*np.random.random((self.num_particles, self.num_dimensions)) - 1.0)
            
    def get_particles(self):
        return np.asarray(self.particles)

    def set_parameters(self, np.ndarray[double, ndim=1, mode="c"] parameters):
        cdef unsigned int i

        if parameters.size != self.num_parameters:
            raise Exception(
                    "Parameters array must be equal to the number "
                    + "parameters")

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

    def __cinit__(self, unsigned int num_particles,
            unsigned int num_dimensions, unsigned int num_parameters,
            double spread=1.0):

        self.wavefunction = new SimpleGaussian(
                num_particles, num_dimensions, num_parameters,
                &self.parameters[0], &self.particles[0, 0])

    def __dealloc__(self):
        pass
        #del self.wavefunction

cdef class PyHamiltonian:
    cdef Hamiltonian *hamiltonian

    def compute_local_energy(self, PyWavefunction wavefunction):
        return self.hamiltonian.compute_local_energy(
                wavefunction.wavefunction)

    def compute_potential_energy(self, PyWavefunction wavefunction):
        return self.hamiltonian.compute_potential_energy(
                wavefunction.wavefunction)

cdef class PyHarmonicOscillator(PyHamiltonian):
    cdef double mass
    cdef double omega

    def __cinit__(self, double mass, double omega):
        self.mass = mass
        self.omega = omega

        self.hamiltonian = new HarmonicOscillator(mass, omega)

    def __dealloc__(self):
        pass
        #del self.hamiltonian

cdef class PyMetropolisAlgorithm:
    cdef MetropolisAlgorithm *method

    def __cinit__(self, unsigned int num_particles):
        self.method = new MetropolisAlgorithm(num_particles)

    def step(self, PyWavefunction wavefunction, double step_length):
        return self.method.step(wavefunction.wavefunction, step_length)

    def run(self, PyWavefunction wavefunction, PyHamiltonian hamiltonian,
            double step_length, unsigned int num_samples):

        return self.method.run(
                wavefunction.wavefunction, hamiltonian.hamiltonian,
                step_length, num_samples)
