#include "wavefunction.h"
#include "math_macros.h"

Wavefunction::Wavefunction(
        unsigned int num_particles,
        unsigned int num_dimensions,
        unsigned int num_parameters,
        double mass,
        double omega,
        double *parameters,
        double *particles)
{
    unsigned int i;

    m_num_particles = num_particles;
    m_num_dimensions = num_dimensions;
    m_num_parameters = num_parameters;

    m_parameters = parameters;

    m_mass = mass;
    m_omega = omega;

    m_particles = new double *[m_num_particles];

    for (i = 0; i < m_num_particles; i++) {
        m_particles[i] = &particles[i*m_num_dimensions];
    }
}

Wavefunction::~Wavefunction()
{
    delete[] m_particles;
}

double Wavefunction::compute_position_squared_sum()
{
    unsigned int i, j;
    double position_squared_sum;

    /* Initialize sum */
    position_squared_sum = 0.0;

    /* Compute the squared sum */
    for (i = 0; i < m_num_particles; i++) {
        for (j = 0; j < m_num_dimensions; j++) {
            position_squared_sum += SQUARE(m_particles[i][j]);
        }
    }

    return position_squared_sum;
}

void Wavefunction::move_particle(double step, unsigned int p_i, unsigned int i)
{
    m_particles[p_i][i] += step;
}
