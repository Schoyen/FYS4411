#pragma once

#include "wavefunction.h"
#include "hamiltonian.h"
#include "monte_carlo_method.h"

class ImportanceMetropolis : public MonteCarloMethod
{
    public:
        using MonteCarloMethod::MonteCarloMethod;

        double greensFraction(
                Wavefunction *wavefunction, double *new_pos, double *old_pos,
                double time_step, double D);

        void initialize();
        bool step(Wavefunction *wavefunction, double step_length);

        //double run(
        //        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
        //        double step_length, unsigned int num_samples);
};
