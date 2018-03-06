#include <valarray>
#include <cmath>

#include "wavefunction.h"
#include "interacting_elliptical_gaussian.h"
#include "math_macros.h"
#include "constants.h"

InteractingEllipticalGaussian::InteractingEllipticalGaussian(
        unsigned int num_particles,
        unsigned int num_dimensions,
        unsigned int num_parameters,
        double mass,
        double omega,
        double beta,
        double *parameters,
        double *particles) :
    Wavefunction(
            num_particles,
            num_dimensions,
            num_parameters,
            mass,
            omega,
            parameters,
            particles)
{
    m_beta = beta;
}

double InteractingEllipticalGaussian::evaluate()
{
    double product;
    unsigned int i, j;

    product = 1.0;

    for (i = 0; i < m_num_particles; i++) {
        product *= evaluate_single_particle_function(i);
        for (j = i + 1; j < m_num_particles; j++) {
            product *= evaluate_correlation_wavefunction(i, j);
        }
    }

    return product;
}

double inline InteractingEllipticalGaussian::evaluate_single_particle_function(
        unsigned int p_i)
{
    double alpha, position_sum, *position;
    unsigned int i;

    alpha = m_parameters[0];

    position = &m_particles[p_i*m_num_dimensions];

    position_sum = 0;

    for (i = 0; i < m_num_dimensions; i++) {
        position_sum +=
            (i != 2) ? SQUARE(position[i]) : SQUARE(m_beta*position[i]);
    }

    return exp(-(alpha*HBAR/(m_mass*m_omega))*position_sum);
}

double InteractingEllipticalGaussian::evaluate_correlation_wavefunction(
        unsigned int p_i, unsigned int p_j)
{
    double a, distance;

    a = m_hard_sphere_radius;

    distance = get_distance_between_particles(p_i, p_j);

    if (distance <= a) {
        return 0.0;
    }

    return 1.0 - a/distance;
}

double InteractingEllipticalGaussian::compute_laplacian()
{
    return 0.0;
}

double InteractingEllipticalGaussian::compute_interaction()
{
    return 0.0;
}

double InteractingEllipticalGaussian::compute_drift_force_component(
        double coordinate)
{
    return 0.0;
}

std::valarray<double>
InteractingEllipticalGaussian::compute_laplacian_variational_gradient()
{
    std::valarray<double> derivative(m_num_parameters);

    derivative = 0;

    return derivative;
}
