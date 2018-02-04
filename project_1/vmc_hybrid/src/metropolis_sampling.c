#include "metropolis_sampling.h"




double perform_metropolis_step(
        wavefunction_t *wavefunction, double step_length)
{
    unsigned int particle_index, dimension_index, i;
    double step, weight,
           old_position[wavefunction->dimensionality];
    particle_t *particle;

    /* Draw a random particle */
    particle_index = arc4random_uniform(wavefunction->num_particles);
    /* Choose random dimension */
    dimension_index = arc4random_uniform(wavefunction->dimensionality);

    /* Create a pointer to the drawn particle */
    particle = &wavefunction->particles[particle_index];

    /* Store the previous position */
    for (i = 0; i < wavefunction->dimensionality; i++) {
        old_position[i] = particle->position[i];
    }

    /* Do a step from [-step_length, step_length) */
    step = step_length*(2.0*RANDOM_UNIFORM_DOUBLE - 1.0);
    /* Propose a new position */
    particle->position[dimension_index] += step;

    /* Compute the ratio between the new and the old position */
    weight = ratio(wavefunction);

    /* Check if we should accept the new state */
    if (weight >= RANDOM_UNIFORM_DOUBLE) {
        /* Store accepted state as the last evaluated value */
        wavefunction->last_value = evaluate_wavefunction(wavefunction);
    } else {
        /* Reset position of particle as we rejected the new state */
        particle->position[dimension_index] = old_position[dimension_index];
    }

    /* Return the local energy */
    return local_energy(wavefunction);
}



double metropolis_sampling(
        wavefunction_t *wavefunction, double step_length,
        unsigned int num_samples)
{
    double energy;
    unsigned int i;

    /* Set initial energy */
    energy = 0;

    /* Perform num_samples metropolis steps */
    for (i = 0; i < num_samples; i++) {
        energy += perform_metropolis_step(wavefunction, step_length);
    }

    /* Return the total energy (without normalization) */
    return energy;
}
