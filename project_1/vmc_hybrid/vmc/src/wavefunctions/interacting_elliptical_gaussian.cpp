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
    unsigned int p_i, p_j;

    product = 1.0;

    for (p_i = 0; p_i < m_num_particles; p_i++) {
        product *= evaluate_single_particle_function(p_i);
        for (p_j = p_i + 1; p_j < m_num_particles; p_j++) {
            product *= evaluate_correlation_wavefunction(p_i, p_j);
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

std::valarray<double>
InteractingEllipticalGaussian::compute_gradient_single_particle_function(
        unsigned int p_k)
{
    std::valarray<double> gradient(m_num_dimensions);
    double alpha, *position;
    unsigned int i;

    alpha = m_parameters[0];

    position = &m_particles[p_k*m_num_dimensions];

    for (i = 0; i < m_num_dimensions; i++) {
        gradient[i] = (i != 2) ? position[i] : (m_beta*position[i]);
    }

    return -2*alpha*gradient;
}

double
InteractingEllipticalGaussian::compute_laplacian_single_particle_function(
        unsigned int p_k)
{
    double alpha, *position, position_sum, laplacian;
    unsigned int i;

    alpha = m_parameters[0];

    position = &m_particles[p_k*m_num_dimensions];

    position_sum = 0;

    for (i = 0; i < m_num_dimensions; i++) {
        position_sum +=
            (i != 2) ? SQUARE(position[i]) : SQUARE(m_beta*position[i]);
    }

    laplacian = -2*alpha*(m_num_dimensions - 1 + m_beta);
    laplacian += 4*SQUARE(alpha)*position_sum;

    return laplacian;
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

std::valarray<double>
InteractingEllipticalGaussian::compute_gradient_correlation_wavefunction(
        unsigned int p_k)
{
    std::valarray<double> gradient(m_num_dimensions);
    double a, *r_k, *r_m, r_km;
    unsigned int p_m, i;

    a = m_hard_sphere_radius;

    r_k = &m_particles[p_k*m_num_dimensions];

    gradient = 0;

    for (p_m = 0; p_m < m_num_particles; p_m++) {
        if (p_m == p_k) {
            continue;
        }

        r_m = &m_particles[p_m*m_num_dimensions];
        r_km = get_distance_between_particles(p_k, p_m);

        for (i = 0; i < m_num_dimensions; i++) {
            gradient[i] += (r_k[i] - r_m[i])*a/(SQUARE(r_km)*(r_km - a));
        }
    }

    return gradient;
}

double
InteractingEllipticalGaussian::compute_laplacian_correlation_wavefunction(
        unsigned int p_k)
{
    double a, *r_k, *r_m, r_km, laplacian;
    unsigned int p_m;

    a = m_hard_sphere_radius;

    r_k = &m_particles[p_k*m_num_dimensions];

    laplacian = 0;

    for (p_m = 0; p_m < m_num_particles; p_m++) {
        if (p_m == p_k) {
            continue;
        }

        r_m = &m_particles[p_m*m_num_dimensions];
        r_km = get_distance_between_particles(p_k, p_m);

        laplacian += (m_num_dimensions - 1)*a/(SQUARE(r_km)*(r_km - a));
        laplacian += (SQUARE(a) - 2*a*r_km)/(SQUARE(r_km)*SQUARE((r_km - a)));
    }

    return laplacian;
}

double InteractingEllipticalGaussian::compute_laplacian(unsigned int p_k)
{
    double laplacian;
    std::valarray<double> gradient_1, gradient_2;
    unsigned int i;

    laplacian = compute_laplacian_single_particle_function(p_k);

    gradient_1 = compute_gradient_single_particle_function(p_k);
    gradient_2 = compute_gradient_correlation_wavefunction(p_k);

    for (i = 0; i < m_num_dimensions; i++) {
        laplacian += 2*gradient_1[i]*gradient_2[i];
        laplacian += SQUARE(gradient_2[i]);
    }

    laplacian += compute_laplacian_correlation_wavefunction(p_k);

    return laplacian;
}

double InteractingEllipticalGaussian::compute_laplacian()
{
    unsigned int p_k;
    double laplacian;

    laplacian = 0;

    for (p_k = 0; p_k < m_num_particles; p_k++) {
        laplacian += compute_laplacian(p_k);
    }

    return laplacian;
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
