#pragma once

#include "wavefunction.h"

class Hamiltonian
{
    public:
        virtual double compute_local_energy(Wavefunction *wavefunction) = 0;
};
