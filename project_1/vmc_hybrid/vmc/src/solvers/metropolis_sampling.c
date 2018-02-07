#include "metropolis_sampling.h"
#include "wavefunction.h"
#include "vmc_macros.h"




double perform_metropolis_step(
        wavefunction_t *wavefunction, hamiltonian_t *hamiltonian,
        double step_length)
{
    unsigned int particle_index, dimension_index, j;
    double step, weight, current_wavefunction,
           old_position[wavefunction->dimensionality];

    /* Draw a random particle */
    particle_index = arc4random_uniform(wavefunction->num_particles);
    /* Choose random dimension */
    dimension_index = arc4random_uniform(wavefunction->dimensionality);

    /* Store the previous position */
    for (j = 0; j < wavefunction->dimensionality; j++) {
        old_position[j] = wavefunction->particles[particle_index][j];
    }

    /* Do a step from [-step_length, step_length) */
    step = step_length*(2.0*RANDOM_UNIFORM_DOUBLE - 1.0);
    /* Propose a new position */
    wavefunction->particles[particle_index][dimension_index] += step;

    /* Evaluate the new wavefunction */
    current_wavefunction = evaluate_wavefunction(wavefunction);

    /* Compute the ratio between the new and the old position */
    weight = SQUARE(current_wavefunction)/SQUARE(wavefunction->last_value);

    /* Check if we should accept the new state */
    if (weight >= RANDOM_UNIFORM_DOUBLE) {
        /* Store accepted state as the last evaluated value */
        wavefunction->last_value = current_wavefunction;
    } else {
        /* Reset position of particle as we rejected the new state */
        wavefunction->particles[particle_index][dimension_index] =
            old_position[dimension_index];
    }

    /* Return the local energy */
    return local_energy(wavefunction, hamiltonian);
}


double metropolis_sampling(
        wavefunction_t *wavefunction, hamiltonian_t *hamiltonian,
        double step_length, unsigned int num_samples)
{
    double energy;
    unsigned int i;

    /* Set initial energy */
    energy = 0;

    /* Perform num_samples metropolis steps */
    for (i = 0; i < num_samples; i++) {
        energy += perform_metropolis_step(wavefunction, hamiltonian,
                step_length);
    }

    /* Return the total energy (without normalization) */
    return energy;
}
