#include <cassert>

#include "hamiltonian.h"
#include "harmonic_oscillator.h"
#include "wavefunction.h"
#include "math_macros.h"
#include "constants.h"

HarmonicOscillator::HarmonicOscillator(double mass, double omega)
{
    assert(mass > 0);
    assert(omega > 0);

    m_mass = mass;
    m_omega = omega;
}

double HarmonicOscillator::compute_local_energy(Wavefunction *wavefunction)
{
    double kinetic_energy, potential_energy;

    kinetic_energy = -SQUARE(HBAR)*wavefunction->compute_laplacian()/(2*m_mass);
    potential_energy = compute_potential_energy(wavefunction);

    return kinetic_energy + potential_energy;
}

double HarmonicOscillator::compute_potential_energy(Wavefunction *wavefunction)
{
    double potential_energy;

    potential_energy = 0.5*m_mass*SQUARE(m_omega);
    potential_energy *= wavefunction->compute_position_squared_sum();

    return potential_energy;
}

double HarmonicOscillator::compute_local_energy_gradient(Wavefunction *wavefunction)
{
    double local_energy_gradient;

    local_energy_gradient = -SQUARE(HBAR)*wavefunction->compute_laplacian_alpha_derivative();
    local_energy_gradient /= 2*m_mass;

    return local_energy_gradient;
}
