import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, PyMetropolisAlgorithm, PySteepestDescentMetropolis, PyImportanceMetropolis
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

num_particles = 1
num_dimensions = 1
num_parameters = 1
spread = 2.0
step_length = 0.05
num_samples = 1000000
mass = 1
omega = 1

alphas = np.linspace(0.1, 1.1, 21)

wavefunction = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
wavefunction.set_parameters(np.array([0.5]))
hamiltonian = PyHarmonicOscillator(mass, omega)
#solver = PyMetropolisAlgorithm(num_particles)
solver = PyImportanceMetropolis(num_particles)

solver.run(wavefunction, hamiltonian, step_length, num_samples)

energies = [] 

for alpha in alphas:

    wavefunction.redistribute()
    wavefunction.set_parameters(np.array([alpha]))

    energy = solver.run(wavefunction, hamiltonian, step_length, num_samples)

    energies.append(energy/num_samples)

    print("alpha: {:.3f}, Energy: {:.3f}".format(alpha, energy/num_samples))

plt.plot(alphas, energies)
plt.show()
