import config

from vmc.interface import PySimpleGaussian, PyHarmonicOscillator, PyMetropolisAlgorithm, PySteepestDescent
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

num_particles = 10
num_dimensions = 2
num_parameters = 1
spread = 5.0
step_length = 0.05
num_samples = 1000000
mass = 1
omega = 1

alphas = np.linspace(0.2, 1.1, 22)

wavefunction = PySimpleGaussian(num_particles, num_dimensions, num_parameters, mass, omega)
hamiltonian = PyHarmonicOscillator(mass, omega)
#solver = PyMetropolisAlgorithm(num_particles)

print(wavefunction.compute_position_squared_sum())
print(wavefunction.compute_alpha_derivative())
#solver = PySteepestDescent()

#solver.step()

energies = [] 
'''
for alpha in alphas:

    wavefunction.redistribute()
    wavefunction.set_parameters(np.array([alpha]))

    energy = solver.run(wavefunction, hamiltonian, step_length, num_samples)

    energies.append(energy/num_samples)

    print("alpha: {:.3f}, Energy: {:.3f}".format(alpha, energy/num_samples))

plt.plot(alphas, energies)
plt.show()

'''

