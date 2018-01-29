#include "metropolis_sampling.h"

double perform_metropolis_step(
        particles_t *particles, double step_length,
        double (*local_energy)(double position),
        double (*ratio)(double new_position, double old_position))
{
    unsigned int particle_index, dimension_index;
    double step, old_position, new_position, weight, delta_energy;
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
    old_position = particle->position[dimension_index];
    /* Propose a new position */
    new_position = old_position + step;

    /* Compute the ratio between the new and the old position */
    weight = (*ratio)(new_position, old_position);

    /* See if the new position should be accepted */
    if (weight >= RANDOM_UNIFORM_DOUBLE) {
        /* Compute the change in energy for the new position */
        delta_energy = (*local_energy)(new_position) \
            - (*local_energy)(old_position);
        /* Update the position of the particle */
        particle->position[dimension_index] = new_position;
    }

    /* Return the change in energy */
    return delta_energy;
}

void metropolis_sampling(void)
{
}
