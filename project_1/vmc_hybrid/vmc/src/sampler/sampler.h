#pragma once

#include <valarray>

#include "monte_carlo_method.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"

class Sampler
{
    private:
        unsigned int m_num_steps;
        unsigned int m_num_accepted_steps;

        double m_energy;
        double m_energy_squared;
        double m_variance;

        /* TODO: This should be a vector */
        double m_position_squared_sum;
        double m_position_energy_sum;

        std::valarray<double> m_energy_gradient;

        unsigned int m_num_local_energies;
        double *m_local_energies;

        Wavefunction *m_wavefunction;
        Hamiltonian *m_hamiltonian;
        MonteCarloMethod *m_solver;
    public:
        Sampler(
                Wavefunction *wavefunction,
                Hamiltonian *hamiltonian,
                MonteCarloMethod *solver,
                unsigned int num_local_energies,
                double *local_energies);

        void initialize()
        {
            m_num_steps = 0;
            m_num_accepted_steps = 0;
            m_energy = 0;
            m_energy_squared = 0;
            m_variance = 0;
            m_energy_gradient = 0;
            m_position_squared_sum = 0;
            m_position_energy_sum = 0;

            if (m_num_local_energies > 0) {
                std::fill(
                        m_local_energies,
                        m_local_energies + m_num_local_energies, 0);
            }
        }

        void normalize()
        {
            m_energy /= m_num_steps;
            m_energy_squared /= m_num_steps;
            m_energy_gradient /= m_num_steps;
            m_position_squared_sum /= m_num_steps;
            m_position_energy_sum /= m_num_steps;
        }

        void sample(unsigned int num_samples, double step_length);

        double get_variance()
        {
            return m_variance;
        }

        double get_energy()
        {
            return m_energy;
        }

        double get_energy_squared()
        {
            return m_energy_squared;
        }

        /* TODO: Lol */
        double get_position_squared_sum()
        {
            return m_position_squared_sum;
        }

        double get_position_energy_sum()
        {
            return m_position_energy_sum;
        }

        double get_acceptance_ratio()
        {
            return ((double) m_num_accepted_steps)/((double) m_num_steps);
        }

        std::valarray<double> get_energy_gradient()
        {
            return m_energy_gradient;
        }
};
