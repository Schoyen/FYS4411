#include <cmath.h>

#include "wavefunction.h"

double local_energy(
        parameters_t *parameters, double position[DIMENSIONALITY])
{
    double position_squared, alpha;
    unsigned int i;

    /* Calculate the position squared */
    position_squared = 0.0;

    for (i = 0; i < DIMENSIONALITY; i++) {
        position_squared += position[i]*position[i];
    }

    /* Fetch the variational parameter */
    alpha = parameters->parameters[i];

    /* Return the local energy */
    return alpha + position_squared*(0.5 - 2*alpha*alpha);
}

double ratio(
        parameters_t *parameters, double new_position[DIMENSIONALITY],
        double old_position[DIMENSIONALITY])
{
    double new_position_squared, old_position_squared, alpha;
    unsigned int i;

    /* Calculate the positions squared */
    new_position_squared = 0.0;
    old_position_squared = 0.0;

    for (i = 0; i < DIMENSIONALITY; i++) {
        new_position_squared += new_position[i]*new_position[i];
        old_position_squared += old_position[i]*old_position[i];
    }

    /* Fetch the variational parameter */
    alpha = parameters->parameters[0];

    /* Return the ratio */
    return exp(-alpha*(new_position_squared - old_position_squared));
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
