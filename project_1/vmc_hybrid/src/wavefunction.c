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
