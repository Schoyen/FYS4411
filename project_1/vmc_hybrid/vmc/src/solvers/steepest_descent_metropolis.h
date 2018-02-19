#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "metropolis_algorithm.h"

class SteepestDescentMetropolis : public MetropolisAlgorithm
{   

    private:

        double m_psiAlphaDerivative1 = 0;
        double m_psiAlphaDerivative2 = 0;

    public:

        using MetropolisAlgorithm::MetropolisAlgorithm;
    
        double steepest_descent(
            Wavefunction *wavefunction, Hamiltonian *hamiltonian,
            double step_length, unsigned int num_samples);

        double run(
                Wavefunction *wavefunction, Hamiltonian *hamiltonian,
                double step_length, unsigned int num_samples);
};
