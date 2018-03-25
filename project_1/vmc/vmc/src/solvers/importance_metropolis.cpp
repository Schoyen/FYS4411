#include "importance_metropolis.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"
#include <iostream>
#include <math.h>
#include <random>


ImportanceMetropolis::ImportanceMetropolis(double diffusion_coefficient)
{
    m_diffusion_coefficient = diffusion_coefficient;
}

ImportanceMetropolis::ImportanceMetropolis(
        double diffusion_coefficient, int seed) :
    MonteCarloMethod(seed)
{
    m_diffusion_coefficient = diffusion_coefficient;
}

bool ImportanceMetropolis::step(Wavefunction *wavefunction, double step_length)
{
    unsigned int num_dimensions, num_particles, p_i, i;
    double previous_wavefunction, current_wavefunction, step, first_vector,
           second_vector, greens_numerator, greens_denominator, weight;

    /* Get the number of dimensions and particles */
    num_dimensions = wavefunction->get_num_dimensions();
    num_particles = wavefunction->get_num_particles();

    /* Store the previous wavefunction for the ratio test */
    previous_wavefunction = wavefunction->evaluate();

    /* Create containers for drift force and positions */
    double old_position[num_dimensions];
    double new_position[num_dimensions];
    double drift_force_old[num_dimensions];
    double drift_force_new[num_dimensions];

    /* Draw a random particle */
    p_i = next_int(0, num_particles - 1);

    /* Store the old position of the random particle */
    wavefunction->copy_particle_position(old_position, p_i);

    /* Compute the drift force for particle p_i */
    wavefunction->compute_drift_force(drift_force_old, p_i);

    /* Propose a new step and move the particle */
    for (i = 0; i < num_dimensions; i++) {
        step = m_diffusion_coefficient*drift_force_old[i]*step_length
            + next_gaussian(0, 1)*sqrt(step_length);

        wavefunction->move_particle(step, p_i, i);
    }

    /* Store the updated particle position */
    wavefunction->copy_particle_position(new_position, p_i);

    /* Evaluate new wavefunction */
    current_wavefunction = wavefunction->evaluate();

    /* Compute the new drift force for particle p_i */
    wavefunction->compute_drift_force(drift_force_new, p_i);

    /* Initialize vector accumulators */
    first_vector = 0;
    second_vector = 0;

    /* Compute the exponential terms in the Green's functions */
    for (i = 0; i < num_dimensions; i++) {
        first_vector += SQUARE(old_position[i] - new_position[i]
                - m_diffusion_coefficient*step_length*drift_force_new[i]);
        second_vector += SQUARE(new_position[i] - old_position[i]
                - m_diffusion_coefficient*step_length*drift_force_old[i]);
    }

    /* Compute the Green's functions */
    greens_numerator =
        exp(-first_vector/(4.0*m_diffusion_coefficient*step_length));
    greens_denominator =
        exp(-second_vector/(4.0*m_diffusion_coefficient*step_length));

    /* Compute the ratio test for the Metropolis-Hastings algorithm */
    weight = (greens_numerator/greens_denominator)
        *SQUARE(current_wavefunction)/SQUARE(previous_wavefunction);

    /* Perform the Metropolis test */
    if (weight >= next_uniform()) {
        return true;
    }

    /* Reset position as we did not accept the state */
    wavefunction->reset_particle_position(old_position, p_i);

    return false;
}
