#pragma once

#include <valarray>

#include "wavefunction.h"

class SimpleGaussian : public Wavefunction
{
    public:
        using Wavefunction::Wavefunction;

        double evaluate();
        double compute_laplacian();
        double compute_drift_force_component(double coordinate);
        std::valarray<double> compute_laplacian_variational_gradient();
};
