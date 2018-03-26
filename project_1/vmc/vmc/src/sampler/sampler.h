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

        double m_r_min = 0;
        double m_r_max = 0;
        double m_bin_step = 0;
        unsigned int m_num_bins = 0;
        double **m_bins = 0;

        std::valarray<double> m_wavefunction_variational_gradient;
        std::valarray<double> m_variational_energy_gradient;

        Wavefunction *m_wavefunction;
        Hamiltonian *m_hamiltonian;
        MonteCarloMethod *m_solver;
    public:
        Sampler(
                Wavefunction *wavefunction,
                Hamiltonian *hamiltonian,
                MonteCarloMethod *solver);

        virtual ~Sampler();

        void initialize()
        {
            unsigned int i, j;

            m_num_steps = 0;
            m_num_accepted_steps = 0;
            m_energy = 0;
            m_energy_squared = 0;
            m_variance = 0;
            m_wavefunction_variational_gradient = 0;
            m_variational_energy_gradient = 0;

            if (m_bins != 0) {
                for (i = 0; i < m_num_bins; i++) {
                    for (j = 0; j < m_wavefunction->get_num_dimensions(); j++) {
                        m_bins[i][j] = 0;
                    }
                }
            }
        }

        void normalize()
        {
            unsigned int i, j;

            m_energy /= m_num_steps;
            m_energy_squared /= m_num_steps;
            m_wavefunction_variational_gradient /= m_num_steps;
            m_variational_energy_gradient /= m_num_steps;

            if (m_bins != 0) {
                for (i = 0; i < m_num_bins; i++) {
                    for (j = 0; j < m_wavefunction->get_num_dimensions(); j++) {
                        m_bins[i][j] /=
                            (m_num_steps*m_wavefunction->get_num_particles());
                    }
                }
            }
        }

        void sample(unsigned int num_samples, double step_length,
                double *local_energies);

        void set_one_body_parameters(
                double r_min, double r_max, unsigned int num_bins, double *bins)
        {
            unsigned int i, num_dimensions;

            m_r_min = r_min;
            m_r_max = r_max;
            m_num_bins = num_bins;
            m_bin_step = (m_r_max - m_r_min)/(m_num_bins - 1);
            m_bins = new double *[m_num_bins];

            num_dimensions = m_wavefunction->get_num_dimensions();
            for (i = 0; i < m_num_bins; i++) {
                m_bins[i] = &bins[i*num_dimensions];
            }
        }

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

        double get_acceptance_ratio()
        {
            return ((double) m_num_accepted_steps)/((double) m_num_steps);
        }

        std::valarray<double> get_variational_parameters_gradient()
        {
            return 2*(m_variational_energy_gradient
                    - m_wavefunction_variational_gradient*m_energy)
                / m_wavefunction->get_num_particles();
        }
};
