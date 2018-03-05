#pragma once

#include <valarray>

#include "hamiltonian.h"
#include "wavefunction.h"

class HarmonicOscillator : public Hamiltonian
{
    public:
        using Hamiltonian::Hamiltonian;

        double compute_local_energy(Wavefunction *wavefunction);
        double compute_potential_energy(Wavefunction *wavefunction);
        std::valarray<double> compute_local_energy_gradient(
                Wavefunction *wavefunction);
};
