#include <valarray>

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
        double *local_energies)
{
    m_wavefunction = wavefunction;
    m_hamiltonian = hamiltonian;
    m_solver = solver;
    m_num_local_energies = num_local_energies;
    m_local_energies = local_energies;

    m_wavefunction_variational_gradient =
        std::valarray<double>(m_wavefunction->get_num_parameters());
    m_variational_energy_gradient =
        std::valarray<double>(m_wavefunction->get_num_parameters());
}

void Sampler::sample(unsigned int num_samples, double step_length)
{
    double current_local_energy;
    unsigned int i;
    std::valarray<double> current_variational_gradient;

    /* Initialize all accumulators and parameters */
    initialize();

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

        /* Get the variational gradient of the wavefunction */
        current_variational_gradient =
            m_wavefunction->compute_variational_gradient();

        /* Accumulate expectation values */
        m_wavefunction_variational_gradient += current_variational_gradient;
        m_variational_energy_gradient +=
            current_variational_gradient*current_local_energy;

        /* Check if we should sample the local energies */
        if (m_num_local_energies != 0) {
            m_local_energies[i] = current_local_energy;
        }
    }

    /* Normalize accumulated energy, energy squared and energy gradient */
    normalize();
    /* Compute the variance */
    m_variance = m_energy_squared - SQUARE(m_energy);
}
