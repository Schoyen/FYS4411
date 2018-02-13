#include <math.h>

#include "wavefunction.h"
#include "simple_gaussian.h"
#include "math_macros.h"
#include "constants.h"
#include <iostream>

SimpleGaussian::SimpleGaussian(
        unsigned int num_particles,
        unsigned int num_dimensions,
        unsigned int num_parameters,
        double mass,
        double omega,
        double *parameters,
        double *particles) :
    Wavefunction(
            num_particles,
            num_dimensions,
            num_parameters,
            parameters,
            particles)
{
    m_mass = mass;
    m_omega = omega;
}

// This includes the division of wavefunction.
double SimpleGaussian::compute_alpha_derivative() 
{
    std::cout << "Started computing.." << std::endl;
    double r_sq = 0;
    std::cout << "Created r_sq" << std::endl;
    std::cout << m_mass << std::endl;
    std::cout << "Number of particles = " << m_num_particles << std::endl;
    std::cout << "Number of dimensions = " << m_num_dimensions << std::endl;
    for (int p = 0; p < m_num_particles; p++) {
        for (int d = 0; d < m_num_dimensions; d++) {
            std::cout << "d = " << d << " p = " << p << std::endl;
            double particle_coordinate = m_particles[p + d*m_num_dimensions];
            r_sq += SQUARE(particle_coordinate);
        }
    }
    return -r_sq;
}

double SimpleGaussian::compute_laplacian()
{
    double alpha, position_squared_sum, laplacian;

    /* Fetch the variational parameter alpha */
    alpha = m_parameters[0];

    /* Compute the position squared sum */
    position_squared_sum = compute_position_squared_sum();

    /* Compute the Laplacian */
    laplacian = (double) -2 * m_num_dimensions * m_num_particles * alpha
        + 4 * SQUARE(alpha) * position_squared_sum;

    return laplacian;
}

double SimpleGaussian::evaluate()
{
    double alpha, evaluated_wavefunction, position_squared_sum;

    /* Fetch the variational parameter alpha */
    alpha = m_parameters[0];

    /* Get the squared sum of the positions */
    position_squared_sum = compute_position_squared_sum();

    /* Compute the wavefunction */
    evaluated_wavefunction =
        exp(-alpha*m_mass*m_omega*position_squared_sum/HBAR);

    /* Return the evaluated result */
    return evaluated_wavefunction;
}
