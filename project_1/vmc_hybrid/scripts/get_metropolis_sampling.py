import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from metropolis_sampling import perform_metropolis, Wavefunction, Particles, \
    normalize_energies

num_particles = 1
num_samples = 1000000
step_length = 0.1
spread = 1.0
num_variational_values = 21
alpha_start = 0.1
alpha_end = 1.4

wavefunction = Wavefunction()
particles = Particles(num_particles)

particles.initialize_particles(spread)

alphas = np.linspace(alpha_start, alpha_end, num_variational_values)

wavefunction.set_parameter_values([0.5])
energies = np.zeros(1)
#energies = np.zeros(alphas.shape)
energies[0] = perform_metropolis(wavefunction, particles, step_length, num_samples)

print (energies[0]/num_particles)
print (normalize_energies(energies, num_samples, num_particles))

#
#for i in tqdm(range(len(alphas))):
#    wavefunction.set_parameter_values([alphas[i]])
#    energies[i] = \
#        perform_metropolis(wavefunction, particles, step_length, num_samples)
#
#energies = normalize_energies(energies, num_samples, num_particles)
#plt.plot(alphas, energies)
#plt.show()
