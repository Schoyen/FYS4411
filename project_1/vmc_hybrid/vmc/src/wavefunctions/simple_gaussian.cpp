#include <math.h>
#include <valarray>

#include "wavefunction.h"
#include "simple_gaussian.h"
#include "math_macros.h"
#include "constants.h"


SimpleGaussian::SimpleGaussian(
        unsigned int num_particles,
        unsigned int num_dimensions,
        double mass,
        double omega,
        double *parameters,
        double *particles) :
    Wavefunction(
            num_particles,
            num_dimensions,
            1, // num_parameters
            mass,
            omega,
            parameters,
            particles)
{
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

std::valarray<double> SimpleGaussian::compute_variational_gradient()
{
    std::valarray<double> gradient(m_num_parameters);

    gradient = compute_position_squared_sum();

    return gradient;
}
