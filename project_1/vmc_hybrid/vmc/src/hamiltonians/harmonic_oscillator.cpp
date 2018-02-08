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

    kinetic_energy = SQUARE(HBAR)*wavefunction->compute_laplacian()/(2*m_mass);
    potential_energy = compute_potential_energy(wavefunction);

    return kinetic_energy + potential_energy;
}

double inline HarmonicOscillator::compute_potential_energy(
        Wavefunction *wavefunction)
{
    double potential_energy;

    potential_energy = 0.5*m_mass*SQUARE(m_omega);
    potential_energy += wavefunction->position_sum_squared();

    return potential_energy;
}

//double HarmonicOscillator::compute_local_energy(Wavefunction *wavefunction)
//{
//    double alpha, position_squared_sum, kinetic_energy, potential_energy;
//    unsigned int i, j;
//
//    /* Fetch the variational parameter */
//    alpha = wavefunction->parameters[0];
//
//    /* Initialize the sum */
//    position_squared_sum = 0.0;
//
//    /* Add all positions */
//    for (i = 0; i < wavefunction->num_particles; i++) {
//        for (j = 0; j < wavefunction->dimensionality; j++) {
//            position_squared_sum += SQUARE(wavefunction->particles[i][j]);
//        }
//    }
//
//    kinetic_energy = -SQUARE(hamiltonian->hbar)/(2.0*hamiltonian->mass) \
//        * (-2*alpha*wavefunction->dimensionality*wavefunction->num_particles \
//            + 4*SQUARE(alpha)*position_squared_sum);
//
//    potential_energy = 0.5*hamiltonian->mass*SQUARE(hamiltonian->omega) \
//        * position_squared_sum;
//
//    /* Return the local energy */
//    return kinetic_energy + potential_energy;
//}
