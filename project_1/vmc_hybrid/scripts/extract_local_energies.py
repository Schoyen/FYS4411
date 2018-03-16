'''
Extract local energies for blocking statistical analysis
'''

import pickle
import numpy as np

from vmc.interface import PyHarmonicOscillator, PyImportanceMetropolis,\
    PySimpleGaussian, PySampler

# Natural units
mass = 1.0
omega = 1.0
hbar = 1.0

# Sampling configuration
step_length = 0.5
alpha = np.array([0.5])
num_particles = 1
num_dimensions = 1
num_parameters = 1

# Monte Carlo classes
solver = PyImportanceMetropolis(num_particles, 0.5)
wavefunction = PySimpleGaussian(num_particles, num_dimensions,\
        mass, omega, spread=step_length)
hamiltonian = PyHarmonicOscillator()
sampler = PySampler(wavefunction, hamiltonian, solver)

number_of_samples = 2000000

sampler.sample(number_of_samples, step_length,\
        num_thermalization_steps=number_of_samples*0.15,\
        sample_local_energies=True)

local_energies = sampler.get_local_energies() / number_of_samples

pickle.dump(local_energies, open("local_energies.p", "wb"))