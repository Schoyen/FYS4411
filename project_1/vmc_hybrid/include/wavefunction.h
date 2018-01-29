#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "wavefunction_config.h"

typedef struct parameters {
    unsigned int num_parameters;
    double parameters[];
} parameters_t;


parameters_t *get_variational_parameters(void);
void free_parameters_struct(parameters_t *parameters);

double local_energy(
        parameters_t *parameters, double position[DIMENSIONALITY]);

double ratio(
        parameters_t *parameters, double new_position[DIMENSIONALITY],
        double old_position[DIMENSIONALITY]);

double wavefunction(
        parameters_t *parameters, double position[DIMENSIONALITY]);

#endif
