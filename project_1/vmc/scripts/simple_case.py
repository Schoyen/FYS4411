import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, \
        PyMetropolisAlgorithm

import numpy as np
import matplotlib.pyplot as plt

import tqdm


num_particles = 50#[1, 10, 100, 500]
num_dimensions = 3
num_parameters = 1
spread = 1.0

step_length = 0.05
num_samples = 1000000

mass = 1.0
omega = 1.0

alpha_min = 0.2
alpha_max = 0.7
num_alphas = 21

alphas = np.linspace(alpha_min, alpha_max, num_alphas).reshape(num_alphas, 1)

wavefunction = PySimpleGaussian(
        num_particles, num_dimensions, num_parameters, spread=spread)
hamiltonian = PyHarmonicOscillator(mass, omega)

method = PyMetropolisAlgorithm(num_particles)

energies = np.zeros(num_alphas)
local_energies_list = []

for i in range(num_alphas):
    wavefunction.set_parameters(alphas[i])
    energies[i] = method.run(
            wavefunction, hamiltonian, step_length, num_samples)

    print ("alphas = {0}\tenergy = {1}".format(
        alphas[i], energies[i]/num_samples))

    wavefunction.redistribute()

plt.plot(alphas.ravel(), energies/num_samples)
plt.show()
