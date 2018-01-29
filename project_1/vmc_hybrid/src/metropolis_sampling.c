#include "metropolis_sampling.h"


double perform_metropolis_step(
        particles_t *particles,
        parameters_t *parameters,
        double step_length)
{
    unsigned int particle_index, dimension_index, i;
    double step, weight, delta_energy;
    double old_position[DIMENSIONALITY], new_position[DIMENSIONALITY];
    particle_t *particle;

    delta_energy = 0;

    /* Draw a random particle */
    particle_index = arc4random_uniform(particles->num_particles);
    /* Choose random dimension */
    dimension_index = arc4random_uniform(DIMENSIONALITY);

    /* Create a pointer to the drawn particle */
    particle = &particles->particles[particle_index];

    /* Do a step from [-step_length, step_length) */
    step = step_length*(2.0*RANDOM_UNIFORM_DOUBLE - 1.0);

    /* Get the previous position */
    for (i = 0; i < DIMENSIONALITY; i++) {
        old_position[i] = particle->position[i];
        new_position[i] = particle->position[i];
    }

    /* Propose a new position */
    new_position[dimension_index] += step;

    /* Compute the ratio between the new and the old position */
    weight = ratio(parameters, new_position, old_position);

    /* See if the new position should be accepted */
    if (weight >= RANDOM_UNIFORM_DOUBLE) {
        /* Compute the change in energy for the new position */
        delta_energy = local_energy(parameters, new_position)
            - local_energy(parameters, old_position);

        /* Update the position of the particle */
        for (i = 0; i < DIMENSIONALITY; i++) {
            particle->position[i] = new_position[i];
        }
    }

    /* Return the change in energy */
    return delta_energy;
}

double metropolis_sampling(
        particles_t *particles,
        parameters_t *parameters,
        double step_length,
        unsigned int num_samples)
{
    double energy;
    unsigned int i;

    /* Compute initial energy */
    energy = 0;

    for (i = 0; i < particles->num_particles; i++) {
        energy += local_energy(parameters, particles->particles[i].position);
    }

    /* Perform num_samples metropolis steps */
    for (i = 0; i < num_samples; i++) {
        energy += perform_metropolis_step(particles, parameters, step_length);
    }

    return energy;
}
