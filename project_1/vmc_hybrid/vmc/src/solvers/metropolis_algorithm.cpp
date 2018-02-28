#include "metropolis_algorithm.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"
#include <iostream>


void MetropolisAlgorithm::initialize()
{
    /* This method does not need an initialization */
}

bool MetropolisAlgorithm::step(Wavefunction *wavefunction, double step_length)
{
    unsigned int particle_index, i, num_dimensions;
    double step, weight, current_wavefunction, previous_wavefunction;

    /* Get evaluated wavefunction prior to moving */
    previous_wavefunction = wavefunction->evaluate();

    /* Get the number of dimensions */
    num_dimensions = wavefunction->get_num_dimensions();

    /* Create a temporary storage for the old position */
    double old_position[num_dimensions];

    /* Draw a random particle */
    particle_index = m_random_particle(m_engine);

    /* Store the previous position */
    wavefunction->copy_particle_position(old_position, particle_index);

    /* Propose a new position */
    for (i = 0; i < num_dimensions; i++) {
        step = step_length*(2.0*m_random_step(m_engine) - 1.0);
        wavefunction->move_particle(step, particle_index, i);
    }

    /* Evaluate the new wavefunction */
    current_wavefunction = wavefunction->evaluate();

    /* Compute the ratio between the new and the old position */
    weight = SQUARE(current_wavefunction)/SQUARE(previous_wavefunction);

    /* Check if we should accept the new state */
    if (weight >= m_random_step(m_engine)) {
        return true;
    } else {
        /* Reset the particle position as we did not accept the state */
        wavefunction->reset_particle_position(old_position, particle_index);
    }

    return false;
}

//double MetropolisAlgorithm::run(
//        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
//        double step_length, unsigned int num_samples)
//{
//    double energy, local_energy;
//    unsigned int i, num_accepted_states;
//
//    /* Set initial number of accepted states */
//    num_accepted_states = 0;
//
//    /* Set initial energy */
//    energy = 0;
//
//    /* Compute initial local energy */
//    local_energy = hamiltonian->compute_local_energy(wavefunction);
//
//    /* Perform num_samples metropolis steps */
//    for (i = 0; i < num_samples; i++) {
//
//        /* Do a step and check if it got accepted */
//        if (step(wavefunction, step_length)) {
//            /* Compute new local energy */
//            local_energy = hamiltonian->compute_local_energy(wavefunction);
//            num_accepted_states++;
//        }
//
//        /* Add local energy */
//        energy += local_energy;
//
//    }
//
//    /* Return the total energy (without normalization) */
//    return energy;
//}
//
//double MetropolisAlgorithm::run_variance(
//        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
//        double step_length, unsigned int num_samples, double *variance)
//{
//    double energy, local_energy, energy_squared;
//    unsigned int i, num_accepted_states;
//
//    /* Set initial number of accepted states */
//    num_accepted_states = 0;
//
//    /* Set initial energy */
//    energy = 0;
//    energy_squared = 0;
//
//    /* Compute initial local energy */
//    local_energy = hamiltonian->compute_local_energy(wavefunction);
//
//    /* Perform num_samples metropolis steps */
//    for (i = 0; i < num_samples; i++) {
//
//        /* Do a step and check if it got accepted */
//        if (step(wavefunction, step_length)) {
//            /* Compute new local energy */
//            local_energy = hamiltonian->compute_local_energy(wavefunction);
//            num_accepted_states++;
//        }
//
//        /* Add local energy */
//        energy += local_energy;
//        energy_squared += SQUARE(local_energy);
//    }
//
//    *variance = energy_squared/num_samples - SQUARE((energy/num_samples));
//
//    /* Return the total energy (without normalization) */
//    return energy;
//}
//
//double MetropolisAlgorithm::run(
//        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
//        double step_length, unsigned int num_samples, double *local_energies)
//{
//    double energy, local_energy;
//    unsigned int i, num_accepted_states;
//
//    /* Set initial number of accepted states */
//    num_accepted_states = 0;
//
//    /* Set initial energy */
//    energy = 0;
//
//    /* Compute initial local energy */
//    local_energy = hamiltonian->compute_local_energy(wavefunction);
//
//    /* Perform num_samples metropolis steps */
//    for (i = 0; i < num_samples; i++) {
//
//        /* Do a step and check if it got accepted */
//        if (step(wavefunction, step_length)) {
//            /* Compute new local energy */
//            local_energy = hamiltonian->compute_local_energy(wavefunction);
//            num_accepted_states++;
//        }
//
//        /* Add local energy */
//        local_energies[i] = local_energy;
//        energy += local_energy;
//    }
//
//    /* Return the total energy (without normalization) */
//    return energy;
//}
