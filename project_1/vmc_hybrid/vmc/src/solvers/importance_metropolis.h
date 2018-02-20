#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "monte_carlo_method.h"

class ImportanceMetropolis : public MonteCarloMethod
{
    public:
        using MonteCarloMethod::MonteCarloMethod;

        bool step(Wavefunction *wavefunction, double step_length);

        double run(
                Wavefunction *wavefunction, Hamiltonian *hamiltonian,
                double step_length, unsigned int num_samples);
};