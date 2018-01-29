#include <cmath.h>

#include "wavefunction.h"

double local_energy(
        parameters_t *parameters, double position[DIMENSIONALITY])
{
    return 0.0;
}

double ratio(
        parameters_t *parameters, double new_position[DIMENSIONALITY],
        double old_position[DIMENSIONALITY])
{
    return 0.0;
}

double wavefunction(
        parameters_t *parameters, double position[DIMENSIONALITY])
{
    double position_squared, alpha;
    unsigned int i;

    /* Calculate the radius/position squared */
    position_squared = 0.0;

    for (i = 0; i < DIMENSIONALITY; i++) {
        position_squared += position[i]*position[i];
    }

    /* Fetch the variational parameter */
    alpha = parameters->parameters[0];

    /* Return the trial wavefunction */
    return exp(-alpha*position_squared);
}
