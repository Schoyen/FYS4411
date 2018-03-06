#include <cmath>

#include "hamiltonian.h"
#include "wavefunction.h"



double Hamiltonian::compute_local_energy(Wavefunction *wavefunction)
{
    double kinetic_energy, potential_energy;

    kinetic_energy = compute_kinetic_energy(wavefunction);
    potential_energy = compute_potential_energy(wavefunction);

    return kinetic_energy + potential_energy;
}

double Hamiltonian::compute_interaction(Wavefunction *wavefunction)
{
    double a, distance;
    const double TOL = 1e-10;
    unsigned int i, j;

    a = wavefunction->get_hard_sphere_radius();

    /* Check if there is a hard sphere */
    if (fabs(a - 0.0) < TOL) {
        return 0.0;
    }

    /* Check all particle combinations */
    for (i = 0; i < wavefunction->get_num_particles(); i++) {
        for (j = i + 1; j < wavefunction->get_num_particles(); j++) {
            distance = wavefunction->get_distance_between_particles(i, j);

            /* Check if the particle hard spheres are overlapping */
            if (fabs(a - distance) <= TOL) {
                /* Abort the loop by returning infinity if any particles are
                 * overlapping */
                return std::numeric_limits<double>::infinity();
            }
        }
    }

    /* No particles are overlapping, return 0 */
    return 0.0;
}
