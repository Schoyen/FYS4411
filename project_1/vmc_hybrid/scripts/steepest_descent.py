'''
    STEEPEST DESCENT (aka Gradient Descent)
    The idea is that a function F would decrease fastest in
    the direction of - Del F(alpha),
    (alpha is a variational parameter).
    We thereby have the iterative scheme:

    alpha[k + 1] = alpha[k] - gamma * Del F(alpha[k])

    if gamma is small enough we have $F(alpha[k]) >= F(alpha[k+1])
    for all k, and hopefulle we have convergense to minimum.
'''

import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, PyImportanceMetropolis, PySampler
from matplotlib import pyplot as plt
import numpy as np

num_particles = 10
num_dimensions = 1
num_parameters = 1
spread = 2.0
step_length = 0.05
num_samples = 1000000

mass = 1
omega = 1

num_local_energies = 100
stride_local_energies = num_samples // num_local_energies


hamiltonian = PyHarmonicOscillator(mass, omega)
wavefunction = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
solver = PyImportanceMetropolis(num_particles)
sampler = PySampler(wavefunction, hamiltonian, solver, num_local_energies, stride_local_energies)

alphas = np.linspace(0.1, 1.1, 11)

'''
for alpha in alphas:
    print("Alpha = ", alpha)
    wavefunction.set_parameters(np.array([alpha]))
    print("Start: ", sampler.get_energy_gradient())
    sampler.sample(num_samples, step_length)
    print("End: ", sampler.get_energy_gradient())
'''

alpha = 0.7 # start

MAX_ITER = 100
iterations = 0
gamma = 0.01
wavefunction.set_parameters(np.array([alpha]))
sampler.sample(num_samples, step_length)
gradient = sampler.get_energy_gradient()
energy = sampler.get_energy()

alphas = np.zeros(MAX_ITER)

while (iterations < MAX_ITER):

    print("Alpha = ", alpha)
    print("Gradient = ", gradient)
    alpha = alpha - gamma * gradient
    wavefunction.set_parameters(np.array([alpha]))
    energy_prev = energy
    energy = sampler.get_energy()
    
    if (energy > energy_prev):
        gamma = gamma*0.95

    sampler.sample(num_samples, step_length)
    gradient = sampler.get_energy_gradient()

    alphas[iterations] = alpha

    iterations += 1

plt.plot(alphas)
plt.show()