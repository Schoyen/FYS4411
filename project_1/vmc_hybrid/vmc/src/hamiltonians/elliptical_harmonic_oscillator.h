#pragma once

#include <valarray>

#include "wavefunction.h"

class EllipticalHarmonicOscillator
{
    private:
        double m_lambda;

    public:
        EllipticalHarmonicOscillator(double lambda);

        double compute_local_energy(Wavefunction *wavefunction);
        double compute_potential_energy(Wavefunction *wavefunction);
        std::valarray<double> compute_local_energy_gradient(
                Wavefunction *wavefunction);
}
