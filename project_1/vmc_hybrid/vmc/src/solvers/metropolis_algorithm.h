#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "monte_carlo_method.h"

class MetropolisAlgorithm : public MonteCarloMethod
{
    public:
        using MonteCarloMethod::MonteCarloMethod;

        void initialize();
        bool step(Wavefunction *wavefunction, double step_length);
};
