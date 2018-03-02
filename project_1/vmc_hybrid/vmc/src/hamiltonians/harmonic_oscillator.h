#pragma once

#include "hamiltonian.h"
#include "wavefunction.h"

class HarmonicOscillator : public Hamiltonian
{
    public:
        HarmonicOscillator(double mass, double omega);
        double compute_local_energy(Wavefunction *wavefunction);
        double compute_potential_energy(Wavefunction *wavefunction);
        double compute_local_energy_gradient(Wavefunction *wavefunction);

    private:
        double m_omega;
        double m_mass;
};
