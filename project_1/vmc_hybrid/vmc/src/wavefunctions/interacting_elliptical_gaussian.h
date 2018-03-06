#pragma once

#include <valarray>

#include "wavefunction.h"

class InteractingEllipticalGaussian : public Wavefunction
{
    private:
        double m_beta;

        double evaluate_single_particle_function(unsigned int p_i);
        std::valarray<double> compute_gradient_single_particle_function(
                unsigned int p_k);
        double compute_laplacian_single_particle_function(unsigned int p_k);

        double evaluate_correlation_wavefunction(
                unsigned int p_i, unsigned int p_j);
        double compute_laplacian(unsigned int p_k);

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
        double compute_drift_force_component(double coordinate);
        std::valarray<double> compute_laplacian_variational_gradient();
};
