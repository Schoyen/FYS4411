#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "metropolis_algorithm.h"

class SteepestDescentMetropolis : public MetropolisAlgorithm
{
    public:

        using MetropolisAlgorithm::MetropolisAlgorithm;
    
        double steepest_descent(
            Wavefunction *wavefunction, Hamiltonian *hamiltonian,
            double step_length, unsigned int num_samples);
};
