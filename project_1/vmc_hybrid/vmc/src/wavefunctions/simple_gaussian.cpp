#include <math.h>
#include <valarray>

#include "wavefunction.h"
#include "simple_gaussian.h"
#include "math_macros.h"
#include "constants.h"



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

std::valarray<double> SimpleGaussian::compute_laplacian_variational_gradient()
{
    double alpha, position_squared_sum;
    std::valarray<double> laplacian_alpha_derivative(m_num_parameters);

    /* Fetch variational parameter alpha */
    alpha = m_parameters[0];

    /* Compute position squared sum */
    position_squared_sum = compute_position_squared_sum();

    /* Compute the alpha derivative of the laplacian */
    laplacian_alpha_derivative[0] =
        (double) -2 * m_num_dimensions * m_num_particles
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
