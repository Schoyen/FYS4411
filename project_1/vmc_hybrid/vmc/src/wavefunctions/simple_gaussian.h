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

        double evaluate();
        double compute_laplacian();
        double compute_drift_force_component(double coordinate);
        std::valarray<double> compute_laplacian_variational_gradient();
};
