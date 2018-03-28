#include <cmath>

#include "double_well.h"
#include "wavefunction.h"
#include "math_macros.h"

DoubleWell::DoubleWell(double radius, unsigned int axis)
{
    m_radius = radius;
    m_axis = axis;
}

double DoubleWell::compute_potential_energy(Wavefunction *wavefunction)
{
    double potential_energy, mass, omega, x_position;
    unsigned int num_particles, num_dimensions, p_i;

    mass = wavefunction->get_mass();
    omega = wavefunction->get_frequency();
    num_particles = wavefunction->get_num_particles();
    num_dimensions = wavefunction->get_num_dimensions();

    potential_energy =
        HarmonicOscillator::compute_potential_energy(wavefunction);

    double position[num_dimensions];

    x_position = 0;
    for (p_i = 0; p_i < num_particles; p_i++) {
        wavefunction->copy_particle_position(position, p_i);
        x_position += fabs(position[m_axis]);
    }

    potential_energy += 0.5*mass*SQUARE(omega)
        *(0.25*SQUARE(m_radius)*num_particles - m_radius*x_position);

    return potential_energy;
}
