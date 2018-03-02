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

double SimpleGaussian::compute_laplacian_alpha_derivative() 
{
    double alpha, position_squared_sum, laplacian_alpha_derivative;

    /* Fetch variational parameter alpha */
    alpha = m_parameters[0];

    /* Compute position squared sum */
    position_squared_sum = compute_position_squared_sum();

    /* Compute the alpha derivative of the laplacian */
    laplacian_alpha_derivative = (double) -2 * m_num_dimensions * m_num_particles
        + 8 * alpha * position_squared_sum;

    return laplacian_alpha_derivative;
}

double SimpleGaussian::compute_drift_force_component(double coordinate)
{
    double alpha; 

    /* Fetch variational parameter alpha */
    alpha = m_parameters[0];

    // Drift force is F = -4*alpha*r_vec
    // return -4*alpha*m_particles[dimension, particle*m_num_dimensions];
    return - 4*alpha*coordinate;
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
