#include <math.h>

#include "wavefunction.h"
#include "simple_gaussian.h"
#include "math_macros.h"

SimpleGaussian::SimpleGaussian(
        unsigned int num_particles,
        unsigned int num_dimensions,
        unsigned int num_parameters,
        double *parameters,
        double *particles) :
    Wavefunction(
            num_particles, num_dimensions, num_parameters, parameters,
            particles)
{
    m_valid_last_value = false;
    evaluate();
}

double SimpleGaussian::compute_laplacian()
{
    double alpha, position_squared_sum;

    /* Fetch the variational parameter alpha */
    alpha = m_parameters[0];

    /* Compute the position squared sum */
    position_squared_sum = compute_position_squared_sum();

    return -2*m_num_dimensions*m_num_particles*alpha
        + 4*SQUARE(alpha)*position_squared_sum;
}

double SimpleGaussian::evaluate()
{
    double alpha, evaluated_wavefunction, position_squared_sum;

    /* Check if we have stored a valid evaluated result and that the particles
     * haven't been moved */
    if (m_valid_last_value && m_valid_position_squared_sum) {
        return m_last_value;
    }

    /* Fetch the variational parameter alpha */
    alpha = m_parameters[0];

    /* Get the squared sum of the positions */
    position_squared_sum = compute_position_squared_sum();

    /* Compute the wavefunction */
    evaluated_wavefunction = exp(-alpha*position_squared_sum);

    /* Update the stored evaluated wavefunction */
    m_last_value = evaluated_wavefunction;
    /* Update the validity of the evaluated wavefunction */
    m_valid_last_value = true;

    /* Return the evaluated result */
    return m_last_value;
}
