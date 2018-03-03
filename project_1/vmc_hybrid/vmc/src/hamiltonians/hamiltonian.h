#pragma once

#include <valarray>

#include "wavefunction.h"

class Hamiltonian
{
    protected:
        double m_a = 0.0043;

    public:
        virtual double compute_local_energy(Wavefunction *wavefunction) = 0;
        virtual double compute_potential_energy(Wavefunction *wavefunction) = 0;
        virtual std::valarray<double> compute_local_energy_gradient(
                Wavefunction *wavefunction) = 0;

        double compute_interaction(Wavefunction *wavefunction);
};
