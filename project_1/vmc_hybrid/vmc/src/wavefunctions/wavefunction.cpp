#include "wavefunction.h"
#include "math_macros.h"

Wavefunction::Wavefunction(
        unsigned int num_particles,
        unsigned int num_dimensions,
        unsigned int num_parameters,
        double *parameters,
        double *particles)
{
    m_num_particles = num_particles;
    m_num_dimensions = num_dimensions;
    m_num_parameters = num_parameters;

    m_parameters = parameters;
    m_particles = particles;

    m_valid_position_squared_sum = false;
    compute_position_squared_sum();
}

double Wavefunction::compute_position_squared_sum()
{
    unsigned int i, j;
    double position_squared_sum;

    /* Check if the last computed value is valid */
    if (m_valid_position_squared_sum) {
        return m_last_position_squared_sum;
    }

    /* Initialize sum */
    position_squared_sum = 0.0;

    /* Compute the squared sum */
    for (i = 0; i < m_num_particles; i++) {
        for (j = 0; j < m_num_dimensions; j++) {
            position_squared_sum += SQUARE(m_particles[j + i*m_num_dimensions]);
        }
    }

    /* Update the stored position sum squared */
    m_last_position_squared_sum = position_squared_sum;
    /* Update the validity of the squared sum */
    m_valid_position_squared_sum = true;

    /* Return the valid squared sum */
    return m_last_position_squared_sum;
}

void Wavefunction::add_step(
        double step, unsigned int particle_index,
        unsigned int coordinate)
{
    m_particles[coordinate + particle_index*m_num_dimensions] += step;
    m_valid_position_squared_sum = false;
}