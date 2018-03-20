/*
 * Hamiltonian with Elliptic Harmonic Oscillator trap potential
 * 
 */

#include <valarray>

#include "wavefunction.h"
#include "elliptical_harmonic_oscillator.h"
#include "hamiltonian.h"
#include "interacting_elliptical_gaussian.h"
#include "math_macros.h"
#include "constants.h"

// Constructor: add the elliptical variational parameter lambda
EllipticalHarmonicOscillator::EllipticalHarmonicOscillator(double lambda)
{
    m_lambda = lambda;
}

double EllipticalHarmonicOscillator::compute_kinetic_energy(
        Wavefunction *wavefunction)
{
    double kinetic_energy;

    kinetic_energy = wavefunction->compute_laplacian();
    kinetic_energy *= -HBAR*wavefunction->get_frequency()/2.0;

    return kinetic_energy;
}

double EllipticalHarmonicOscillator::compute_potential_energy(
        Wavefunction *wavefunction)
{
    double potential_energy, position[wavefunction->get_num_dimensions()],
        position_sum;
    unsigned int i, p_k;

    position_sum = 0;

    // Iterating over all particles. 
    for (p_k = 0; p_k < wavefunction->get_num_particles(); p_k++) {
        wavefunction->copy_particle_position(position, p_k);

        // Iterating over every particle coordinate (dimensions)
        for (i = 0; i < wavefunction->get_num_dimensions(); i++) {
            
            // For the third dimension - add the elliptical perturbation
            position_sum +=
                (i != 2) ? SQUARE(position[i]) : SQUARE(position[i]*m_lambda);
        }
    }

    potential_energy = HBAR*wavefunction->get_frequency()/2.0;
    potential_energy *= position_sum;

    return potential_energy;
}
