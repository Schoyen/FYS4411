#include <math.h>
#include <stdio.h>

#include "wavefunction.h"
#include "vmc_macros.h"

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
