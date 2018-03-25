#pragma once

#include <valarray>

#include "hamiltonian.h"
#include "wavefunction.h"
#include "interacting_elliptical_gaussian.h"

class EllipticalHarmonicOscillator : public Hamiltonian
{
    private:
        double m_lambda;

    public:
        EllipticalHarmonicOscillator(double lambda);

        double compute_kinetic_energy(Wavefunction *wavefunction);
        double compute_potential_energy(Wavefunction *wavefunction);
};
