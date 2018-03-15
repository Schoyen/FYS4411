#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "monte_carlo_method.h"

class ImportanceMetropolis : public MonteCarloMethod
{
    private:
        double m_diffusion_coefficient = 0.5;
        double m_time_step;

    public:
        ImportanceMetropolis(double time_step, double diffusion_coefficient);
        ImportanceMetropolis(
                double time_step, double diffusion_coefficient, int seed);

        bool step(Wavefunction *wavefunction, double step_length);
};
