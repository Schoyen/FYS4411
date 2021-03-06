/*
 * Gaussian wavefunction for use with elliptic traph
 *
 */

#include <valarray>
#include <cmath>

#include "wavefunction.h"
#include "interacting_elliptical_gaussian.h"
#include "math_macros.h"
#include "constants.h"

InteractingEllipticalGaussian::InteractingEllipticalGaussian(
        unsigned int num_particles,
        unsigned int num_dimensions,
        double mass,
        double omega,
        double beta,
        double radius,
        double *parameters,
        double *particles) :
    Wavefunction(
            num_particles,
            num_dimensions,
            1, //num_parameters
            mass,
            omega,
            parameters,
            particles)
{
    m_beta = beta;
    m_hard_sphere_radius = radius;
}

double InteractingEllipticalGaussian::evaluate()
{
    double product;
    unsigned int p_i, p_j;

    product = 1.0;

    // Wavefunction evaluation, single-particle and interaction contributions
    for (p_i = 0; p_i < m_num_particles; p_i++) {
        product *= evaluate_single_particle_function(p_i);
        for (p_j = p_i + 1; p_j < m_num_particles; p_j++) {
            product *= evaluate_correlation_wavefunction(p_i, p_j);
        }
    }

    return product;
}

// Private method. Single-particle contributions.
double inline InteractingEllipticalGaussian::evaluate_single_particle_function(
        unsigned int p_i)
{
    double alpha, position_sum, *position;
    unsigned int i;

    alpha = m_parameters[0];

    position = m_particles[p_i];

    position_sum = 0;

    for (i = 0; i < m_num_dimensions; i++) {

        // The z-coordinates will be perturbed.
        position_sum +=
            (i != 2) ? SQUARE(position[i]) : (m_beta*SQUARE(position[i]));
    }

    return exp(-alpha*position_sum);
}

// Private method. Interaction (correlation) contributions.
double InteractingEllipticalGaussian::evaluate_correlation_wavefunction(
        unsigned int p_i, unsigned int p_j)
{
    double a, distance;

    a = m_hard_sphere_radius;

    distance = get_distance_between_particles(p_i, p_j);

    // Removing "illegal" position contributions
    if (distance <= a) {
        return 0.0;
    }

    return 1.0 - a/distance;
}

// Method to compute laplacian for all particles.
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

// PRIVAET method to compute laplacian for a single particle.
double InteractingEllipticalGaussian::compute_laplacian(unsigned int p_k)
{
    double laplacian;
    std::valarray<double> gradient_spf, gradient_cw;
    unsigned int i;

    // Call to helper functions
    laplacian = compute_laplacian_single_particle_function(p_k);
    gradient_spf = compute_gradient_single_particle_function(p_k);
    gradient_cw = compute_gradient_correlation_wavefunction(p_k);

    // Adding up contributions for the particular in question.
    for (i = 0; i < m_num_dimensions; i++) {
        laplacian += 2*gradient_spf[i]*gradient_cw[i];
        laplacian += SQUARE(gradient_cw[i]);
    }

    laplacian += compute_laplacian_correlation_wavefunction(p_k);

    return laplacian;
}

// Private method. Helper. Compute SPF gradient.
std::valarray<double>
InteractingEllipticalGaussian::compute_gradient_single_particle_function(
        unsigned int p_k)
{
    std::valarray<double> gradient(m_num_dimensions);
    double alpha, *position;
    unsigned int i;

    alpha = m_parameters[0];

    position = m_particles[p_k];

    // Conditional z-coordinate perturbation.
    for (i = 0; i < m_num_dimensions; i++) {
        gradient[i] = (i != 2) ? position[i] : (m_beta*position[i]);
    }

    return -2*alpha*gradient;
}

// Private method. Helper. Compute SPF laplacian.
double
InteractingEllipticalGaussian::compute_laplacian_single_particle_function(
        unsigned int p_k)
{
    double alpha, *position, position_sum, laplacian;
    unsigned int i;

    alpha = m_parameters[0];

    position = m_particles[p_k];

    position_sum = 0;

    // Conditional z-coordinate perturbation.
    for (i = 0; i < m_num_dimensions; i++) {
        position_sum +=
            (i != 2) ? SQUARE(position[i]) : SQUARE(m_beta*position[i]);
    }

    laplacian = -2*alpha*(m_num_dimensions - 1 + m_beta);
    laplacian += 4*SQUARE(alpha)*position_sum;

    return laplacian;
}

// Private method. Computes CW gradient for scpecific particle.
std::valarray<double>
InteractingEllipticalGaussian::compute_gradient_correlation_wavefunction(
        unsigned int p_k)
{
    std::valarray<double> gradient(m_num_dimensions);
    double a, *r_k, *r_m, r_km;
    unsigned int p_m, i;

    a = m_hard_sphere_radius;

    r_k = m_particles[p_k];

    gradient = 0;

    // Compute all particle contributions. Ignore if self (p_m == p_k).
    for (p_m = 0; p_m < m_num_particles; p_m++) {
        if (p_m == p_k) {
            continue;
        }

        // Getting distance b/w particles.j
        r_m = m_particles[p_m];
        r_km = get_distance_between_particles(p_k, p_m);

        // Adding contribution.
        for (i = 0; i < m_num_dimensions; i++) {
            gradient[i] += (r_k[i] - r_m[i])*a/(SQUARE(r_km)*(r_km - a));
        }
    }

    return gradient;
}

// Private method. Computes CW laplacian for specific particle.
double
InteractingEllipticalGaussian::compute_laplacian_correlation_wavefunction(
        unsigned int p_k)
{
    double a, r_km, laplacian;
    unsigned int p_m;

    a = m_hard_sphere_radius;

    laplacian = 0;

    // Compute all particle contributions. Ignore if self (p_m == p_k).
    for (p_m = 0; p_m < m_num_particles; p_m++) {
        if (p_m == p_k) {
            continue;
        }

        // Get distance b/w particles.
        r_km = get_distance_between_particles(p_k, p_m);

        // Adding contribution.
        laplacian += (m_num_dimensions - 1)*a/(SQUARE(r_km)*(r_km - a));
        laplacian += (SQUARE(a) - 2*a*r_km)/(SQUARE(r_km)*SQUARE((r_km - a)));
    }

    return laplacian;
}

// Public method. Full gradient for particular particle.
void InteractingEllipticalGaussian::compute_gradient(
        double *gradient, unsigned int p_i)
{
    unsigned int i;
    std::valarray<double> spf_gradient, corr_gradient;

    // Calling private helper methods.
    spf_gradient = compute_gradient_single_particle_function(p_i);
    corr_gradient = compute_gradient_correlation_wavefunction(p_i);

    // Contribution from all particles to particle at p_i
    for (i = 0; i < m_num_dimensions; i++) {
        gradient[i] = spf_gradient[i] + corr_gradient[i];
    }
}

// Public method. Full variational gradient for particular particle.
std::valarray<double>
InteractingEllipticalGaussian::compute_variational_gradient()
{
    unsigned int p_i, i;
    std::valarray<double> gradient(m_num_parameters);
    double particle_gradient;

    particle_gradient = 0;

    // Iterating over all particles and coordinates. Adding up contributions.
    for (p_i = 0; p_i < m_num_particles; p_i++) {
        for (i = 0; i < m_num_dimensions; i++) {
            particle_gradient +=
                (i != 2) ? SQUARE(m_particles[p_i][i])
                : (m_beta*SQUARE(m_particles[p_i][i]));
        }
    }

    gradient = -particle_gradient;

    return gradient;
}
