import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, PyMetropolisAlgorithm, PySteepestDescentMetropolis, PyImportanceMetropolis
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

num_particles = 100
num_dimensions = 3
num_parameters = 1
spread = 5.0
step_length = 0.05
num_samples = 1000000
mass = 1
omega = 1

no_of_alphas = 11
alphas = np.linspace(0.3, 0.7, no_of_alphas)

wavefunction_1 = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
wavefunction_2 = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
hamiltonian_1 = PyHarmonicOscillator(mass, omega)
hamiltonian_2 = PyHarmonicOscillator(mass, omega)
regular_metropolis_solver = PyMetropolisAlgorithm(num_particles)
importance_metropolis_solver = PyImportanceMetropolis(num_particles)

energies = np.zeros((2, no_of_alphas)) 

for i in range(no_of_alphas):

    wavefunction_1.redistribute()
    wavefunction_2.redistribute()
    wavefunction_1.set_parameters(np.array([alphas[i]]))
    wavefunction_2.set_parameters(np.array([alphas[i]]))

    regular_energy    = regular_metropolis_solver.run(wavefunction_1, hamiltonian_1, step_length, num_samples)
    importance_energy = importance_metropolis_solver.run(wavefunction_2, hamiltonian_2, step_length, num_samples) 

    regular_energy = regular_energy / num_samples
    importance_energy = importance_energy / num_samples

    energies[0, i] = regular_energy 
    energies[1, i] = importance_energy 

    print("alpha: {:.3f}, Reg. metr. energy: {:.3f}, Imp. samp. energy. {:.3f}.".format(alphas[i], regular_energy, importance_energy))


plt.plot(alphas, energies[0, :], '--r', label="Regular metropolis")
plt.plot(alphas, energies[1, :], '--b', label="With importance sampling")
plt.title(str(num_dimensions) + "D, N = " + str(num_particles))
plt.ylabel("Energy")
plt.xlabel("alpha")
plt.legend()
plt.show()
