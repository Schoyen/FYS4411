import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, PyMetropolisAlgorithm
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

num_particles = 10
num_dimensions = 3
num_parameters = 1
spread = 1.0
step_length = 0.05
num_samples = 1000
mass = 1
omega = 1

alphas = np.linspace(0, 2, 11)

wavefunction = PySimpleGaussian(num_particles, num_dimensions, num_parameters)
hamiltonian = PyHarmonicOscillator(mass, omega)
solver = PyMetropolisAlgorithm(num_particles)

for alpha in alphas:

    wavefunction.set_parameters(np.array([alpha]))

    solver.run(wavefunction, hamiltonian, step_length, num_samples)

    print("alpha: ", alpha, "Wavefunction value: ", wavefunction.evaluate())





