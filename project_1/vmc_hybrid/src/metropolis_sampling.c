#include "metropolis_sampling.h"

double perform_metropolis_step(particles_t *particles, double step_length)
{
    unsigned int particle_index, dimension_index;
    double step, proposed_position;
    particle_t *particle;

    /* Draw a random particle */
    particle_index = arc4random_uniform(particles->num_particles);
    /* Choose random dimension */
    dimension_index = arc4random_uniform(DIMENSIONALITY);

    /* Create a pointer to the drawn particle */
    particle = &particles->particles[particle_index];

    /* Do a step from [-step_length, step_length) */
    step = step_length*(2.0*RANDOM_UNIFORM_DOUBLE - 1.0);
    /* Propose a new position */
    proposed_position = particle->position[dimension_index] + step;

    return step;
}

void metropolis_sampling(void)
{
}
