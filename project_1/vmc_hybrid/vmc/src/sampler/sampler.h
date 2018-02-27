#pragma once

#include "monte_carlo_method.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"

class Sampler
{
    private:
        unsigned int m_num_steps;
        unsigned int m_accepted_steps;

        double m_energy;
        double m_energy_squared;
        double m_variance;
        double *m_local_energies;

        Wavefunction *m_wavefunction;
        Hamiltonian *m_hamiltonian;
    public:
        Sampler(
                Wavefunction *wavefunction,
                Hamiltonian *hamiltonian,
                double *local_energies);

        void initialize()
        {
            m_num_steps = 0;
            m_accepted_steps = 0;
            m_energy = 0;
            m_energy_squared = 0;
            m_variance = 0;
        }

        void sample(unsigned int num_samples, double step_length);

        void compute_variance()
        {
            m_variance =
                m_energy_squared/num_samples - SQUARE((m_energy/m_num_steps));
        }

        double get_variance()
        {
            return m_variance;
        }

        double get_energy()
        {
            return m_energy;
        }
};
