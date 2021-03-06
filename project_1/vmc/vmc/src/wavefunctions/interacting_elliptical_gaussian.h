#pragma once

#include <valarray>

#include "wavefunction.h"

class InteractingEllipticalGaussian : public Wavefunction
{

   // Private helper functions
    private:
        double m_beta;

        double evaluate_single_particle_function(unsigned int p_i);
        std::valarray<double> compute_gradient_single_particle_function(
                unsigned int p_k);
        double compute_laplacian_single_particle_function(unsigned int p_k);

        double evaluate_correlation_wavefunction(
                unsigned int p_i, unsigned int p_j);
        std::valarray<double> compute_gradient_correlation_wavefunction(
                unsigned int p_k);
        double compute_laplacian_correlation_wavefunction(unsigned int p_k);

        double compute_laplacian(unsigned int p_k);

    public:
        InteractingEllipticalGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                double mass,
                double omega,
                double beta,
                double radius,
                double *parameters,
                double *particles);

        double get_beta()
        {
            return m_beta;
        }

        double evaluate();
        double compute_laplacian();
        void compute_gradient(double *gradient, unsigned int p_i);
        std::valarray<double> compute_variational_gradient();
};
