#pragma once

#include "wavefunction.h"

class SimpleGaussian : public Wavefunction
{
    protected:
        double m_mass;
        double m_omega;
        //unsigned int m_num_particles;
        //unsigned int m_num_dimensions;
    public:
        SimpleGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double mass,
                double omega,
                double *parameters,
                double *particles);

        double evaluate();
        double compute_alpha_derivative();
        double compute_laplacian();
};
