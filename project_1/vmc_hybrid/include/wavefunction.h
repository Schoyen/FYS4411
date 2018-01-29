#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "particles.h"

typedef struct parameters {
    unsigned int num_parameters;
    double parameters[];
} parameters_t;


double local_energy(
        parameters_t *parameters, double position[DIMENSIONALITY]);

double ratio(
        parameters_t *parameters, double new_position[DIMENSIONALITY],
        double old_position[DIMENSIONALITY]);

double wavefunction(
        parameters_t *parameters, double position[DIMENSIONALITY]);

#endif
