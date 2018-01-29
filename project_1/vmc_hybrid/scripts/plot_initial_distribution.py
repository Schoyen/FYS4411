from particles import Particles
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

system = Particles(50)
system.initialize_particles(3.0)

particle_array = system.get_particles()

ax = plt.figure().add_subplot(111, projection="3d")
ax.scatter(particle_array[:, 0], particle_array[:, 1], particle_array[:, 2])
plt.show()
