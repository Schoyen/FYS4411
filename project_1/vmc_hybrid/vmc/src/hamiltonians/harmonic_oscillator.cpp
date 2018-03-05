#include <cassert>
#include <valarray>

#include "hamiltonian.h"
#include "harmonic_oscillator.h"
#include "wavefunction.h"
#include "math_macros.h"
#include "constants.h"

double HarmonicOscillator::compute_local_energy(Wavefunction *wavefunction)
{
    double kinetic_energy, potential_energy;

    kinetic_energy = wavefunction->compute_laplacian();
    kinetic_energy *= -SQUARE(HBAR)/(2*wavefunction->get_mass());

    potential_energy = compute_potential_energy(wavefunction);

    return kinetic_energy + potential_energy;
}

double HarmonicOscillator::compute_potential_energy(Wavefunction *wavefunction)
{
    double potential_energy;

    potential_energy = wavefunction->compute_position_squared_sum();
    potential_energy *= 0.5*wavefunction->get_mass();
    potential_energy *= SQUARE((wavefunction->get_frequency()));

    return potential_energy;
}

std::valarray<double> HarmonicOscillator::compute_local_energy_gradient(
        Wavefunction *wavefunction)
{
    std::valarray<double> local_energy_gradient;

    local_energy_gradient =
        wavefunction->compute_laplacian_variational_gradient();
    local_energy_gradient *= -SQUARE(HBAR)/(2*wavefunction->get_mass());

    return local_energy_gradient;
}
