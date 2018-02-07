#include <stdlib.h>
#include <stdio.h>

#include "wavefunction.h"


void allocate_variational_parameters(wavefunction_t *wavefunction)
{
    /* Allocate memory for parameters */
    wavefunction->parameters =
        (double *) calloc(wavefunction->num_parameters, sizeof(double));

    /* Check if the allocation succeeded */
    if (!wavefunction->parameters) {
        fprintf(stderr, "Allocation of parameters failed\n");
        exit(EXIT_FAILURE);
    }
}

void free_variational_parameters(wavefunction_t *wavefunction)
{
    /* Free parameter array memory */
    free(wavefunction->parameters);
}

void allocate_particles(wavefunction_t *wavefunction)
{
    unsigned int i;

    /* Allocate memory for each particle */
    wavefunction->particles =
        (double **) calloc(wavefunction->num_particles, sizeof(double *));

    /* Check if the allocation was successful */
    if (!wavefunction->particles) {
        fprintf(stderr, "Allocation of particle rows failed\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate entire position matrix as one long contiguous array */
    wavefunction->particles[0] =
        (double *) calloc(
                wavefunction->num_particles*wavefunction->dimensionality,
                sizeof(double));

    /* Check if the allocation passed */
    if (!wavefunction->particles[0]) {
        fprintf(stderr, "Allocation of all particle positions failed\n");
        exit(EXIT_FAILURE);
    }

    /* Set pointers for each particle */
    for (i = 1; i < wavefunction->num_particles; i++) {
        wavefunction->particles[i] =
            wavefunction->particles[0] + i*wavefunction->dimensionality;
    }
}

void free_particles(wavefunction_t *wavefunction)
{
    /* Free position array */
    free(wavefunction->particles[0]);
    /* Free particle array */
    free(wavefunction->particles);
}
