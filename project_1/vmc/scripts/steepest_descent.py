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

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, \
    PyImportanceMetropolis, PySampler, PyEllipticalHarmonicOscillator, \
    PyInteractingEllipticalGaussian
from matplotlib import pyplot as plt
import numpy as np

num_particles = 10
num_dimensions = 3
step_length = 1.0
num_samples = int(1e5)
num_thermalization_steps = 0

mass = 1
omega = 1

num_local_energies = 0

beta = _lambda = 2.82843
radius = 0.043


#hamiltonian = PyHarmonicOscillator()
hamiltonian = PyEllipticalHarmonicOscillator(_lambda)
wavefunction = PyInteractingEllipticalGaussian(
        num_particles, num_dimensions, mass, omega, beta, radius,
        spread=step_length)
#wavefunction = PySimpleGaussian(num_particles, num_dimensions, mass, omega)
solver = PyImportanceMetropolis()
sampler = PySampler(wavefunction, hamiltonian, solver)

alpha = np.array([1.0]) # start

MAX_ITER = 100
iterations = 0
gamma = 0.001

alphas = np.zeros(MAX_ITER)

gradient_prev = 1e10

while (iterations < MAX_ITER):

    wavefunction.set_parameters(alpha)
    sampler.sample(
            num_samples, step_length,
            num_thermalization_steps=num_thermalization_steps)

    print("Alpha = ", alpha)
    gradient = sampler.get_parameter_gradient()
    print("Gradient = ", gradient)
    alpha = alpha - gamma * gradient
    wavefunction.redistribute()

    #if (gradient*gradient > gradient_prev*gradient_prev):
    #    print ("Boink")
    #    gamma = gamma * 0.9

    gradient_prev = gradient

    mask = alpha < 0
    alpha[mask] = -alpha[mask]

    alphas[iterations] = alpha[0]
    iterations += 1

plt.plot(alphas)
plt.show()
