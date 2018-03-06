#pragma once

#include <valarray>

#include "wavefunction.h"

class InteractingEllipticalGaussian : public Wavefunction
{
    private:
        double m_beta;

        double evaluate_single_particle_function(unsigned int i);
        double correlation_wavefunction(unsigned int i, unsigned int j);

    public:
        InteractingEllipticalGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double mass,
                double omega,
                double beta,
                double *parameters,
                double *particles);

        double get_beta()
        {
            return m_beta;
        }

        double evaluate();
        double compute_laplacian();
        double compute_interaction();
        double compute_drift_force_component(double coordinate);
        std::valarray<double> compute_laplacian_variational_gradient();
};
