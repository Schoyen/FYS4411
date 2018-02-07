#include "hamiltonian.h"
#include "wavefunction.h"
#include "vmc_macros.h"

double local_energy(wavefunction_t *wavefunction, hamiltonian_t *hamiltonian)
{
    double alpha, position_squared_sum, kinetic_energy, potential_energy;
    unsigned int i, j;

    /* Fetch the variational parameter */
    alpha = wavefunction->parameters[0];

    /* Initialize the sum */
    position_squared_sum = 0.0;

    /* Add all positions */
    for (i = 0; i < wavefunction->num_particles; i++) {
        for (j = 0; j < wavefunction->dimensionality; j++) {
            position_squared_sum += SQUARE(wavefunction->particles[i][j]);
        }
    }

    kinetic_energy = -SQUARE(hamiltonian->hbar)/(2.0*hamiltonian->mass) \
        * (-2*alpha*wavefunction->dimensionality*wavefunction->num_particles \
            + 4*SQUARE(alpha)*position_squared_sum);

    potential_energy = 0.5*hamiltonian->mass*SQUARE(hamiltonian->omega) \
        * position_squared_sum;

    /* Return the local energy */
    return kinetic_energy + potential_energy;
}

