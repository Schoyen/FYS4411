'''
    PROBLEM 1 B: Deveoping the code. 

    Compute the ground state for a Spherical Harmonic Oscillator.

    Use both numerical and analytic computation of local energy. 
    Only parameter is alpha. Compare CPU time difference. 

    Number of particles.  N = {1, 10, 100, 500}
    Number of dimensions. d = {1, 2, 3}

    Compare with exact answer. -> Is the ground state energy correct?
'''

# Importing necessary tools
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np 
import time

from vmc.interface import PyHarmonicOscillator, PyMetropolisAlgorithm, \
    PySimpleGaussian, PySimpleGaussianNumerical, PySampler


# Natural units
mass  = 1.0
omega = 1.0
hbar  = 1.0

# Initial distribution spread
spread      = 1.0

# Pretermined step particles (walkers) can move for each iteration
step_length = 0.05 

# Parameter space, variational parameter alpha (include alpha = 0.5)
num_alphas = 11
alpha = np.linspace(0.3, 0.7, num_alphas).reshape(num_alphas, 1)

# Particle-, dimension- and parameter configuration
num_particles  = [1, 10, 100, 500]
num_dimensions = [1, 2, 3]
num_parameters = 1

# Storage
energies_analytic  = np.zeros((len(num_particles), len(num_dimensions), num_alphas))
energies_numerical = np.zeros((len(num_particles), len(num_dimensions), num_alphas))
cpu_time_analytic  = np.zeros((len(num_particles), len(num_dimensions), num_alphas))
cpu_time_numerical = np.zeros((len(num_particles), len(num_dimensions), num_alphas))

# THE BIG HORRIBLE THING
for i in range(len(num_particles)):
    for j in range(len(num_dimensions)):

        print(" ---------------- ")
        print("Number of particles: {}, {}D".format(num_particles[i], num_dimensions[j]))

        # Monte Carlo specific parameters
        solver = PyMetropolisAlgorithm(num_particles[i])
        analytic_wfn  = PySimpleGaussian(num_particles[i], num_dimensions[j], num_parameters, mass, omega, spread=spread)
        numerical_wfn = PySimpleGaussianNumerical(num_particles[i], num_dimensions[j], num_parameters, mass, omega, spread=spread)
        hamiltonian  = PyHarmonicOscillator()
        analytic_sampler  = PySampler(analytic_wfn, hamiltonian, solver, 0)
        numerical_sampler = PySampler(numerical_wfn,  hamiltonian, solver, 0)
        num_samples = int(1200 * num_particles[i])

        # Iterating over alphas
        for k in range(num_alphas):

            # Setting variational parameter
            analytic_wfn.set_parameters(alpha[k])
            numerical_wfn.set_parameters(alpha[k])

            # Sampling, analytic wavefunction
            analytic_start_time = time.time()
            analytic_sampler.sample(num_samples, step_length, num_thermalization_steps=int(0.15*num_samples))
            analytic_stop_time = time.time()

            # Storing analytic time and energy
            cpu_time_analytic[i, j, k] += analytic_stop_time - analytic_start_time
            energies_analytic[i, j, k] += analytic_sampler.get_energy()
             
            # Sampling, numeric wavefunction
            numerical_start_time = time.time()
            numerical_sampler.sample(num_samples, step_length, num_thermalization_steps=int(0.15*num_samples))
            numerical_stop_time = time.time()

            # Storing numerical time and energy
            cpu_time_numerical[i, j, k] += numerical_stop_time - numerical_start_time
            energies_numerical[i, j, k] += numerical_sampler.get_energy()

            # Redistribute walkers/particles
            analytic_wfn.redistribute()
            numerical_wfn.redistribute()

            print("Alpha = {:5.3f}, Analytic E = {:9.5f}, CPU Time: {:8.5f}, Numerical E = {:9.5f}, CPU Time: {:8.5f}"
                    .format(alpha.ravel()[k], energies_analytic[i, j, k], cpu_time_analytic[i, j, k], energies_numerical[i, j, k], cpu_time_numerical[i, j, k]))