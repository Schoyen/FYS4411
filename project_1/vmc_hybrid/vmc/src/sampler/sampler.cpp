#include "sampler.h"
#include "monte_carlo_method.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"

Sampler::Sampler(
        Wavefunction *wavefunction,
        Hamiltonian *hamiltonian,
        MonteCarloMethod *solver,
        unsigned int num_local_energies,
        unsigned int stride_local_energies,
        double *local_energies)
{
    m_wavefunction = wavefunction;
    m_hamiltonian = hamiltonian;
    m_solver = solver;
    m_num_local_energies = num_local_energies;
    m_stride_local_energies = stride_local_energies;
    m_local_energies = local_energies;
}

void Sampler::sample(unsigned int num_samples, double step_length)
{
    double current_local_energy;
    unsigned int i;

    /* Initialize all accumulators and parameters */
    initialize();

    /* Initialize the solver */
    m_solver->initialize();

    /* Store the number of samples which will be performed */
    m_num_steps = num_samples;
    /* Compute initial local energy */
    current_local_energy = m_hamiltonian->compute_local_energy(m_wavefunction);

    /* Perform m_num_steps metropolis steps */
    for (i = 0; i < m_num_steps; i++) {

        /* Do a step and check if it got accepted */
        if (m_solver->step(m_wavefunction, step_length)) {
            /* Compute new local energy */
            current_local_energy =
                m_hamiltonian->compute_local_energy(m_wavefunction);
            m_num_accepted_steps++;
        }

        /* Add local energy */
        m_energy += current_local_energy;
        /* Add local energy squared */
        m_energy_squared += SQUARE(current_local_energy);
        /* Add local energy gradient */
        /* TODO: Should be a vector */
        m_energy_gradient +=
            m_hamiltonian->compute_local_energy_gradient(m_wavefunction);

        /* Check if we should sample the local energies and if we are at the
         * correct iteration */
        if (m_num_local_energies != 0 && i % m_stride_local_energies) {
            m_local_energies[i / m_stride_local_energies] =
                current_local_energy;
        }
    }

    /* Normalize accumulated energy, energy squared and energy gradient */
    normalize_energies();
    /* Compute the variance */
    m_variance = m_energy_squared - SQUARE(m_energy);
}
