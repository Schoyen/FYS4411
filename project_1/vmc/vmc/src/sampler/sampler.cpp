#include <valarray>

#include "sampler.h"
#include "monte_carlo_method.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"

Sampler::Sampler(
        Wavefunction *wavefunction,
        Hamiltonian *hamiltonian,
        MonteCarloMethod *solver)
{
    m_wavefunction = wavefunction;
    m_hamiltonian = hamiltonian;
    m_solver = solver;

    m_wavefunction_variational_gradient =
        std::valarray<double>(m_wavefunction->get_num_parameters());
    m_variational_energy_gradient =
        std::valarray<double>(m_wavefunction->get_num_parameters());
}

Sampler::~Sampler()
{
    if (m_num_bins != 0) {
        delete[] m_bins;
    }
}

void Sampler::sample(unsigned int num_samples, double step_length,
        double *local_energies)
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
        if (local_energies != 0) {
            local_energies[i] = current_local_energy;
        }

        /* Check if we should sample one body densities */
        if (m_bins != 0) {
            unsigned int p_i, i, bin_index, num_particles, num_dimensions;

            num_dimensions = m_wavefunction->get_num_dimensions();
            num_particles = m_wavefunction->get_num_particles();

            double position[num_dimensions];
            for (p_i = 0; p_i < num_particles; p_i++) {
                m_wavefunction->copy_particle_position(position, p_i);

                for (i = 0; i < num_dimensions; i++) {
                    position[i] = fabs(position[i]);
                    if (position[i] < m_r_min || position[i] > m_r_max) {
                        continue;
                    }

                    bin_index = (int) floor((position[i] - m_r_min)/m_bin_step);
                    m_bins[bin_index][i] += 1;
                }
            }
        }
    }
    /* Normalize accumulated energy, energy squared and energy gradient */
    normalize();

    /* Compute the variance */
    m_variance = m_energy_squared - SQUARE(m_energy);
    m_variance /= m_num_steps;
}
