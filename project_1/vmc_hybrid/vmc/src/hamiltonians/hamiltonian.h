#pragma once

#include "wavefunction.h"

class Hamiltonian
{
    public:
        virtual double compute_local_energy(Wavefunction *wavefunction) = 0;
        virtual double compute_potential_energy(Wavefunction *wavefunction) = 0;
        virtual double compute_local_energy_gradient(Wavefunction *wavefunction) = 0;
};
