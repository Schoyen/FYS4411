#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "wavefunction.h"


/* This system uses a single variational parameter alpha */
void allocate_variational_parameters(wavefunction_t *wavefunction)
{
    /* Allocate memory for alpha */
    wavefunction->parameters = (double *) malloc(sizeof(double));

    /* Check if the allocation succeeded */
    if (!wavefunction->parameters) {
        fprintf(stderr, "Allocation of parameters failed\n");
        exit(EXIT_FAILURE);
    }

    /* Set the number of parameters in the system, i.e., only alpha */
    wavefunction->num_parameters = 1;
    /* Set initial value for alpha */
    wavefunction->parameters[0] = 1.0;
}

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

double local_energy_total(
        parameters_t *parameters, particles_t *particles)
{
    double energy;
    unsigned int i;

    /* Initialize energy */
    energy = 0.0;

    /* Add all the energy contribution */
    for (i = 0; i < particles->num_particles; i++) {
        energy += local_energy(parameters, particles->particles[i].position);
    }

    /* Return the total local energy (unormalized) */
    return energy;
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
    return exp(-2*alpha*(new_position_squared - old_position_squared));
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
