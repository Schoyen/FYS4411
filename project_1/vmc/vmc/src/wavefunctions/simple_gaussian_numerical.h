#pragma once

#include "wavefunction.h"
#include "simple_gaussian.h"

class SimpleGaussianNumerical : public SimpleGaussian
{
    private:
        double m_h;

    public:
        SimpleGaussianNumerical(
                unsigned int num_particles,
                unsigned int num_dimensions,
                double mass,
                double omega,
                double h,
                double *parameters,
                double *particles);

        double compute_laplacian();
};
