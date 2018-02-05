#include <math.h>
#include <stdio.h>

#include "wavefunction.h"


double local_energy(wavefunction_t *wavefunction)
{
    double alpha, position_squared_sum;
    unsigned int i, j;

    /* Fetch the variational parameter */
    alpha = wavefunction->parameters[0];

    /* Initialize the sum */
    position_squared_sum = 0.0;

    // TODO: Store the parameters of the Hamiltonian, e.g., omega, mass

    /* Add all positions */
    for (i = 0; i < wavefunction->num_particles; i++) {
        for (j = 0; j < wavefunction->dimensionality; j++) {
            position_squared_sum += SQUARE(wavefunction->particles[i][j]);
        }
    }

    /* Return the local energy */
    return (-2*SQUARE(alpha) + 0.5)*position_squared_sum
        + wavefunction->dimensionality*wavefunction->num_particles*alpha;
}

double ratio(wavefunction_t *wavefunction)
{
    double current_value;

    /* Get the current value of the wavefunction */
    current_value = evaluate_wavefunction(wavefunction);

    /* Return the ratio of the old and the new wavefunction */
    return SQUARE(current_value)/SQUARE(wavefunction->last_value);
}

double evaluate_wavefunction(wavefunction_t *wavefunction)
{
    double position_squared_sum, alpha;
    unsigned int i, j;

    /* Initialize position squared sum */
    position_squared_sum = 0.0;

    /* Calculate position squared sum */
    for (i = 0; i < wavefunction->num_particles; i++) {
        for (j = 0; j < wavefunction->dimensionality; j++) {
            position_squared_sum += SQUARE(wavefunction->particles[i][j]);
        }
    }

    /* Fetch the variational parameter */
    alpha = wavefunction->parameters[0];

    /* Return the trial wavefunction */
    return exp(-alpha*position_squared_sum);
}
