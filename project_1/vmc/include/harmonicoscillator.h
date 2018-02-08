#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian {
    public:
        HarmonicOscillator(System* system, double omega);
        double computeLocalEnergy(Particles* walkers);

    private:
        double m_omega = 0;
};