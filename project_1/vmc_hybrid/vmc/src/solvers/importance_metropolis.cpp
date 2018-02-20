#include "importance_metropolis.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"
#include <iostream>

bool ImportanceMetropolis::step(Wavefunction *wavefunction, double step_length)
{
    /*
        This method is potentially very similar to the other step methods.
        Uses a different ratio to determine if step is accepted:

        ( G(old_pos, new_pos, delta_t) * | psi_T(new)|^2 ) /
        ( G(new_pos, old_pos, delta_t) * | psi_T(old)|^2 )
    */

    std::cout << "is this thing on?" << std::endl;
    return false;
}

double ImportanceMetropolis::run(
        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
        double step_length, unsigned int num_samples)
{
    double energy, local_energy;
    unsigned int i, num_accepted_states;

    /* Set initial number of accepted states */
    num_accepted_states = 0;

    /* Set initial energy */
    energy = 0;

    /* Compute initial local energy */
    local_energy = hamiltonian->compute_local_energy(wavefunction);

    /* Perform num_samples metropolis steps */
    for (i = 0; i < num_samples; i++) {

        /* Do a step and check if it got accepted */
        if (step(wavefunction, step_length)) {
            /* Compute new local energy */
            local_energy = hamiltonian->compute_local_energy(wavefunction);
            num_accepted_states++;
        }

        /* Add local energy */
        energy += local_energy;

    }

    /* Return the total energy (without normalization) */
    return energy;
}