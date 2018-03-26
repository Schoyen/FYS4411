#include <valarray>

#include "hamiltonian.h"
#include "harmonic_oscillator.h"
#include "wavefunction.h"
#include "math_macros.h"
#include "constants.h"



double HarmonicOscillator::compute_kinetic_energy(Wavefunction *wavefunction)
{
    double kinetic_energy;

    kinetic_energy = wavefunction->compute_laplacian();
    kinetic_energy *= -SQUARE(HBAR)/(2*wavefunction->get_mass());

    return kinetic_energy;
}

double HarmonicOscillator::compute_potential_energy(Wavefunction *wavefunction)
{
    double potential_energy;

    potential_energy = wavefunction->compute_position_squared_sum();
    potential_energy *= 0.5*wavefunction->get_mass();
    potential_energy *= SQUARE((wavefunction->get_frequency()));

    return potential_energy;
}
