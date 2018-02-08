#pragma once

#include "wavefunction.h"

class SimpleGaussian : public Wavefunction
{
    public:
        SimpleGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double *parameters,
                double *particles);

        double evaluate();
        double compute_laplacian();
};
