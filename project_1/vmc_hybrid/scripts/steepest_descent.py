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

num_particles = 20
num_dimensions = 1
num_parameters = 1
spread = 1.0
step_length = 0.01
num_samples = int(1e6)
num_thermalization_steps = int(1e6)

mass = 1
omega = 1

num_local_energies = 100
stride_local_energies = num_samples // num_local_energies


hamiltonian = PyHarmonicOscillator()
wavefunction = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
solver = PyImportanceMetropolis(num_particles)
sampler = PySampler(wavefunction, hamiltonian, solver, num_local_energies, stride_local_energies)

alpha = 1.0
wavefunction.set_parameters(np.array([alpha]))
gradient = 0
experiments = 100
for i in range(1, experiments + 1):
    sampler.sample(num_samples, step_length, num_thermalization_steps)
    gradient += sampler.get_energy_gradient()
    wavefunction.redistribute()
    print ("Intermediate gradient: {0}".format(gradient/i))

print ("Final gradient: {0}".format(gradient/experiments))
__import__('sys').exit()

alphas = np.linspace(0.1, 1.1, 11)

for alpha in alphas:
    print("Alpha = ", alpha)
    wavefunction.set_parameters(np.array([alpha]))
    print("Start: ", sampler.get_energy_gradient())
    sampler.sample(num_samples, step_length)
    print("End: ", sampler.get_energy_gradient())
    wavefunction.redistribute()

alpha = 0.5 # start

MAX_ITER = 100
iterations = 0
iters_since_last_improvement = 0
gamma = 1.0
wavefunction.set_parameters(np.array([alpha]))
sampler.sample(
        num_samples, step_length,
        num_thermalization_steps=num_thermalization_steps)
gradient = sampler.get_energy_gradient()# / num_particles
energy = sampler.get_energy()

alphas = np.zeros(MAX_ITER)

while (iterations < MAX_ITER):

    print("Alpha = ", alpha)
    print("Gradient = ", gradient)
    alpha = alpha - gamma * gradient
    wavefunction.set_parameters(np.array([alpha]))
    wavefunction.redistribute()
    energy_prev = energy
    energy = sampler.get_energy()
    gradient_prev = gradient
    sampler.sample(
            num_samples, step_length,
            num_thermalization_steps=num_thermalization_steps)
    gradient = sampler.get_energy_gradient()# /num_particles

    if (gradient*gradient > gradient_prev*gradient_prev):
        if (iters_since_last_improvement >= 10):
            gamma = 0.01
            print("Stuck!")
        else:
            gamma = gamma * 0.9

        iters_since_last_improvement += 1
    else:
        iters_since_last_improvement = 0

    if alpha < 0:
        alpha = -alpha

    alphas[iterations] = alpha
    iterations += 1

plt.plot(alphas)
plt.show()
