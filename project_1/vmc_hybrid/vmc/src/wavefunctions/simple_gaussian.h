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

        double compute_gradient_component(
                unsigned int p_i, unsigned int i)
        {
            double alpha;

            /* Fetch variational parameter alpha */
            alpha = m_parameters[0];

            return -2*alpha*m_particles[p_i][i];
        }

        double evaluate();
        double compute_laplacian();
        std::valarray<double> compute_laplacian_variational_gradient();
};
