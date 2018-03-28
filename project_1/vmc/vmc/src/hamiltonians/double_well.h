#pragma once

#include "harmonic_oscillator.h"
#include "wavefunction.h"

class DoubleWell : public HarmonicOscillator
{
    private:
        double m_radius = 0;
        unsigned int m_axis = 0;

    public:
        DoubleWell(double radius, unsigned int axis);

        double compute_potential_energy(Wavefunction *wavefunction);
};
