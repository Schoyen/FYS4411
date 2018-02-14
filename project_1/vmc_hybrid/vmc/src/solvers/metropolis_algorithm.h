#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "monte_carlo_method.h"

class MetropolisAlgorithm : public MonteCarloMethod
{

    protected:
        double m_psiAlphaDerivative1 = 0;
        double m_psiAlphaDerivative2 = 0;

    public:
        using MonteCarloMethod::MonteCarloMethod;

        bool step(Wavefunction *wavefunction, double step_length);
        double run(
                Wavefunction *wavefunction, Hamiltonian *hamiltonian,
                double step_length, unsigned int num_samples);
        double run_variance(
                Wavefunction *wavefunction, Hamiltonian *hamiltonian,
                double step_length, unsigned int num_samples,
                double *variance);
        double run(
                Wavefunction *wavefunction, Hamiltonian *hamiltonian,
                double step_length, unsigned int num_samples,
                double *local_energies);
};
