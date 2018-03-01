import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, PyMetropolisAlgorithm, PyImportanceMetropolis, PySampler
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

# !@$#!?...
num_local_energies = 100
stride_local_energies = num_samples // num_local_energies

wavefunction_1 = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
wavefunction_2 = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
hamiltonian_1 = PyHarmonicOscillator(mass, omega)
hamiltonian_2 = PyHarmonicOscillator(mass, omega)
regular_metropolis_solver = PyMetropolisAlgorithm(num_particles)
importance_metropolis_solver = PyImportanceMetropolis(num_particles)
sampler_1 = PySampler(wavefunction_1, hamiltonian_1, regular_metropolis_solver, num_local_energies, stride_local_energies)
sampler_2 = PySampler(wavefunction_2, hamiltonian_2, importance_metropolis_solver, num_local_energies, stride_local_energies)

energies = np.zeros((2, no_of_alphas)) 

for i in range(no_of_alphas):

    wavefunction_1.set_parameters(np.array([alphas[i]]))
    wavefunction_2.set_parameters(np.array([alphas[i]]))

    sampler_1.sample(num_samples, step_length)
    sampler_2.sample(num_samples, step_length)

    regular_energy = sampler_1.get_energy() 
    importance_energy = sampler_2.get_energy()

    energies[0, i] = regular_energy 
    energies[1, i] = importance_energy 

    wavefunction_1.redistribute()
    wavefunction_2.redistribute()

    print("alpha: {:.3f}, Reg. metr. energy: {:.3f}, Imp. samp. energy. {:.3f}.".format(alphas[i], regular_energy, importance_energy))


plt.plot(alphas, energies[0, :], '--r', label="Regular metropolis")
plt.plot(alphas, energies[1, :], '--b', label="With importance sampling")
plt.title(str(num_dimensions) + "D, N = " + str(num_particles))
plt.ylabel("Energy")
plt.xlabel("alpha")
plt.legend()
plt.show()
