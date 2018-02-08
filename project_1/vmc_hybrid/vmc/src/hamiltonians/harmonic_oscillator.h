#pragma once

#include "hamiltonian.h"
#include "wavefunction.h"

class HarmonicOscillator : public Hamiltonian
{
    public:
        HarmonicOscillator(double mass, double omega);
        double compute_local_energy(Wavefunction *wavefunction);

    private:
        double m_omega = 0;
        double m_mass = 1;
};
