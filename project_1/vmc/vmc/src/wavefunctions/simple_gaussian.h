#pragma once

#include <valarray>

#include "wavefunction.h"

class SimpleGaussian : public Wavefunction
{
    public:
        SimpleGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                double mass,
                double omega,
                double *parameters,
                double *particles);

        void compute_gradient(double *gradient, unsigned int p_i)
        {
            double alpha;
            unsigned int i;

            alpha = m_parameters[0];

            for (i = 0; i < m_num_dimensions; i++) {
                gradient[i] = -2*alpha*m_particles[p_i][i];
            }
        }

        double evaluate();
        double compute_laplacian();
        std::valarray<double> compute_variational_gradient();
};
