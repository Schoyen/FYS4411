from metropolis_sampling import Wavefunction, perform_varying_metropolis, \
        normalize_energies
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

num_particles = 1
dimensionality = 3
num_parameters = 1
spread = 1.0
step_length = 0.5
num_samples = 100000

parameters = np.array([[alpha] for alpha in np.linspace(0.1, 2, 101)])

wavefunction = Wavefunction(num_particles, dimensionality, num_parameters)
wavefunction.initialize_wavefunction(spread)

particle_array = wavefunction.get_particles()


# One dimensional plotting for DIMENSIONALITY = 1
#plt.scatter(np.linspace(0, 50, len(particle_array)), particle_array)
#plt.show()

# Three dimensional plotting for DIMENSIONALITY = 3
ax = plt.figure().add_subplot(111, projection="3d")
ax.scatter(particle_array[:, 0], particle_array[:, 1], particle_array[:, 2])
plt.show()


energies = perform_varying_metropolis(
        wavefunction, step_length, num_samples, parameters)

energies = normalize_energies(energies, num_samples, num_particles)

particle_array = wavefunction.get_particles()
ax = plt.figure().add_subplot(111, projection="3d")
ax.scatter(particle_array[:, 0], particle_array[:, 1], particle_array[:, 2])
plt.show()

plt.plot(parameters[:, 0], energies)
plt.show()
