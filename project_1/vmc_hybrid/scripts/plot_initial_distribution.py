from metropolis_sampling import Particles
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

system = Particles(50)
system.initialize_particles(3.0)

particle_array = system.get_particles()

# One dimensional plotting for DIMENSIONALITY = 1
#plt.scatter(np.linspace(0, 50, len(particle_array)), particle_array)
#plt.show()

# Three dimensional plotting for DIMENSIONALITY = 3
ax = plt.figure().add_subplot(111, projection="3d")
ax.scatter(particle_array[:, 0], particle_array[:, 1], particle_array[:, 2])
plt.show()
