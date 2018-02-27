#pragma once

#include "wavefunction.h"

class SimpleGaussian : public Wavefunction
{
    protected:
        double m_mass;
        double m_omega;
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
        double compute_laplacian();
        //double compute_drift_force_component(unsigned int dimension, unsigned int particle);
        double compute_drift_force_component(double coordinate);
};
