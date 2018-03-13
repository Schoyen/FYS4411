#include "importance_metropolis.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"
#include <iostream>
#include <math.h>
#include <random>


ImportanceMetropolis::ImportanceMetropolis(
        double time_step, double diffusion_coefficient)
{
    m_time_step = time_step;
    m_diffusion_coefficient = diffusion_coefficient;
}

ImportanceMetropolis::ImportanceMetropolis(
        double time_step, double diffusion_coefficient, int seed) :
    MonteCarloMethod(seed)
{
    m_time_step = time_step;
    m_diffusion_coefficient = diffusion_coefficient;
}

// Function to compute the Green's functions fraction
double ImportanceMetropolis::greensFraction(
        Wavefunction *wavefunction, double* new_pos, double* old_pos)
{
    double result;
    double drift_force_old, drift_force_new;
    double G_numerator = 0;
    double G_denominator = 0;
    unsigned int i, num_dimensions;

    double first_vector = 0;
    double second_vector = 0;

    num_dimensions = wavefunction->get_num_dimensions();

    for (i = 0; i < num_dimensions; i++) {

        drift_force_old = wavefunction->compute_drift_force_component(old_pos[i]);
        drift_force_new = wavefunction->compute_drift_force_component(new_pos[i]);

        first_vector  += SQUARE(old_pos[i] - new_pos[i]
                - m_diffusion_coefficient*m_time_step*drift_force_new);
        second_vector += SQUARE(new_pos[i] - old_pos[i]
                - m_diffusion_coefficient*m_time_step*drift_force_old);
    }

    G_numerator = exp(-first_vector / (4*m_diffusion_coefficient*m_time_step));
    G_denominator = exp(-second_vector / (4*m_diffusion_coefficient*m_time_step));

    result = G_numerator / G_denominator;

    return result;
}

bool ImportanceMetropolis::step(Wavefunction *wavefunction, double step_length)
{
    /*
        This method is potentially very similar to the other step methods.
        Uses a different ratio to determine if step is accepted:

        ( G(old_pos, new_pos, delta_t) * | psi_T(new)|^2 ) /
        ( G(new_pos, old_pos, delta_t) * | psi_T(old)|^2 )

        And then there was ths whole deal of what G is..
        What is nice: Both G's have a coefficient that is cancelled.
    */

    unsigned int num_dimensions, num_particles, i, particle_index;
    double step, weight, current_wavefunction, previous_wavefunction;

    //// The Diffusion coefficient. I don't know what it is supposed to be.
    //// Hardcoding
    //double D = 0.5;

    //// Is this so?
    //double time_step = step_length;

    // Getting dimensionality
    num_dimensions = wavefunction->get_num_dimensions();

    // Get the number of particles
    num_particles = wavefunction->get_num_particles();

    // Storing the current wavefunction in prev (current will be overwritten)
    previous_wavefunction = wavefunction->evaluate();

    // Temporary storage for old position
    double old_pos[num_dimensions];

    // Temporary storage for new position
    double new_pos[num_dimensions];

    // draw some particle
    particle_index = next_int(0, num_particles - 1);

    // Store old (this) position
    wavefunction->copy_particle_position(old_pos, particle_index);

    /*
        Update position:
        r_new = r_prev + D*F(r)*time_step + N(0,1)*sqrt(time_step)
        -> step = r_new - r_prev = D*F*time_step + N *sqrt(time_step)
        N is normal dist. stoc. var.
        D is some diffusion coefficient. WHAT IS IT!?!?!?1
    */

    double drift_force = 0;

    // Popose new position
    for (i = 0; i < num_dimensions; i++) {
        drift_force = wavefunction->compute_drift_force_component(old_pos[i]);

        step =
            m_diffusion_coefficient*drift_force*m_time_step
            + next_gaussian(0, 1)*sqrt(m_time_step);

        wavefunction->move_particle(step, particle_index, i);
    }

    // Store new position
    wavefunction->copy_particle_position(new_pos, particle_index);

    // Evaluate new wavefunction
    current_wavefunction = wavefunction->evaluate();

    // Compute new weight
    double greens_fraction = greensFraction(wavefunction, new_pos, old_pos);

    // Acceptance weight
    weight = greens_fraction
             * (SQUARE(current_wavefunction) / SQUARE(previous_wavefunction));

    // The same kind of test as before
    if (weight >= next_uniform()) {
        return true;
    } else {
        wavefunction->reset_particle_position(old_pos, particle_index);
    }

    return false;
}
