#include "wavefunction.h"
#include "simple_gaussian.h"
#include "simple_gaussian_numerical.h"
#include "math_macros.h"

SimpleGaussianNumerical::SimpleGaussianNumerical(
        unsigned int num_particles,
        unsigned int num_dimensions,
        unsigned int num_parameters,
        double mass,
        double omega,
        double h,
        double *parameters,
        double *particles) :
    SimpleGaussian(
            num_particles,
            num_dimensions,
            num_parameters,
            mass,
            omega,
            parameters,
            particles)
{
    m_h = h;
}

double SimpleGaussianNumerical::compute_laplacian()
{
    double laplacian, central, forward, backward, h_squared;
    unsigned int i, j;

    h_squared = SQUARE(m_h);
    central = evaluate();
    laplacian = 0;

    for (i = 0; i < m_num_particles; i++) {
        for (j = 0; j < m_num_dimensions; j++) {
            move_particle(m_h, i, j);
            forward = evaluate();

            move_particle(-2*m_h, i, j);
            backward = evaluate();

            move_particle(m_h, i, j);
            laplacian += (forward - 2*central + backward)/h_squared;
        }
    }

    return laplacian;
}
