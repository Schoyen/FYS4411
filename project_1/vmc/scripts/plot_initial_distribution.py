import config

from vmc.interface import PySimpleGaussian
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

num_particles = 1
num_dimensions = 3
num_parameters = 1
spread = 2.0
step_length = 0.5
num_samples = 100000

wavefunction = PySimpleGaussian(num_particles, num_dimensions, num_parameters)

particle_array = wavefunction.get_particles()
ax = plt.figure().add_subplot(111, projection="3d")
ax.scatter(particle_array[:, 0], particle_array[:, 1], particle_array[:, 2])
plt.show()

wavefunction.set_parameters(np.array([1.0]))

print (wavefunction.evaluate())
print (wavefunction.compute_laplacian())

wavefunction.set_parameters(np.array([0.5]))

print (wavefunction.evaluate())
print (wavefunction.compute_laplacian())
