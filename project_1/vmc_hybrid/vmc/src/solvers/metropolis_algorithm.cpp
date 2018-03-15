#include "metropolis_algorithm.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"
#include <iostream>



bool MetropolisAlgorithm::step(Wavefunction *wavefunction, double step_length)
{
    unsigned int particle_index, i, num_dimensions, num_particles;
    double step, weight, current_wavefunction, previous_wavefunction;

    /* Get evaluated wavefunction prior to moving */
    previous_wavefunction = wavefunction->evaluate();

    /* Get the number of dimensions */
    num_dimensions = wavefunction->get_num_dimensions();
    /* Get the number of particles */
    num_particles = wavefunction->get_num_particles();

    /* Create a temporary storage for the old position */
    double old_position[num_dimensions];

    /* Draw a random particle */
    particle_index = next_int(0, num_particles - 1);

    /* Store the previous position */
    wavefunction->copy_particle_position(old_position, particle_index);

    /* Propose a new position */
    for (i = 0; i < num_dimensions; i++) {
        step = step_length*next_uniform(-1.0, 1.0);
        wavefunction->move_particle(step, particle_index, i);
    }

    /* Evaluate the new wavefunction */
    current_wavefunction = wavefunction->evaluate();

    /* Compute the ratio between the new and the old position */
    weight = SQUARE(current_wavefunction)/SQUARE(previous_wavefunction);

    /* Check if we should accept the new state */
    if (weight >= next_uniform()) {
        return true;
    }

    /* Reset the particle position as we did not accept the state */
    wavefunction->reset_particle_position(old_position, particle_index);

    return false;
}
