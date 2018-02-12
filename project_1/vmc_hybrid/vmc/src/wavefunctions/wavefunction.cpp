#include "wavefunction.h"
#include "math_macros.h"
#include "constants.h"

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

double Wavefunction::compute_numerical_laplacian()
{
    double laplacian, central, forward, backward, h_squared;
    unsigned int i, j;

    h_squared = SQUARE(DELTA_H);
    central = evaluate();
    laplacian = 0;

    for (i = 0; i < m_num_particles; i++) {
        for (j = 0; j < m_num_dimensions; j++) {
            move_particle(DELTA_H, i, j);
            forward = evaluate();

            move_particle(-2*DELTA_H, i, j);
            backward = evaluate();

            move_particle(DELTA_H, i, j);
            laplacian += (forward - 2*central + backward)/h_squared;
        }
    }

    return laplacian;
}

void Wavefunction::move_particle(
        double step, unsigned int particle_index,
        unsigned int coordinate)
{
    m_particles[coordinate + particle_index*m_num_dimensions] += step;
    m_valid_position_squared_sum = false;
}
