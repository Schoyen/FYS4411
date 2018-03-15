#pragma once

#include <valarray>

#include "hamiltonian.h"
#include "wavefunction.h"

class HarmonicOscillator : public Hamiltonian
{
    public:
        using Hamiltonian::Hamiltonian;

        double compute_kinetic_energy(Wavefunction *wavefunction);
        double compute_potential_energy(Wavefunction *wavefunction);
};
