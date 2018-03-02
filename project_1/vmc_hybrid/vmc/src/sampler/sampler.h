#pragma once

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
        unsigned int m_num_local_energies;
        unsigned int m_stride_local_energies;
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
                unsigned int stride_local_energies,
                double *local_energies);

        void initialize()
        {
            m_num_steps = 0;
            m_num_accepted_steps = 0;
            m_energy = 0;
            m_energy_squared = 0;
            m_variance = 0;

            if (m_num_local_energies > 0) {
                std::fill(
                        m_local_energies,
                        m_local_energies + m_num_local_energies, 0);
            }
        }

        void normalize_energies()
        {
            m_energy /= m_num_steps;
            m_energy_squared /= m_num_steps;
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

        double get_ratio_of_accepted_steps()
        {
            return ((double) m_num_accepted_steps)/((double) m_num_steps);
        }

        double get_energy_gradient();
};
