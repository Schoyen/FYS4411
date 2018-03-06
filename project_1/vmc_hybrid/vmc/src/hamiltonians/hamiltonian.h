#pragma once

#include <valarray>

#include "wavefunction.h"

class Hamiltonian
{
    public:
        virtual double compute_kinetic_energy(Wavefunction *wavefunction) = 0;
        virtual double compute_potential_energy(Wavefunction *wavefunction) = 0;
        virtual std::valarray<double> compute_local_energy_gradient(
                Wavefunction *wavefunction) = 0;

        double compute_local_energy(Wavefunction *wavefunction);
        double compute_interaction(Wavefunction *wavefunction);
};
