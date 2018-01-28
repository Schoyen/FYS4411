#include "metropolis_sampling.h"

double perform_metropolis_step(particles_t *particles, double step_length)
{
    /* Draw a random particle */
    unsigned int particle_index = arc4random_uniform(particles->num_particles);
    /* Choose random dimension */
    unsigned int dimension_index = arc4random_uniform(DIMENSIONALITY);
    /* Create a pointer to the drawn particle */
    particle_t *particle = &particles->particles[particle_index];
    /* Do a step from [-step_length, step_length) */
    double step = step_length*(2.0*RANDOM_UNIFORM_DOUBLE - 1.0);
    /* Propose a new position */
    double proposed_position = particle->position[dimension_index] + step;

    return step;
}

void metropolis_sampling(void)
{
}
