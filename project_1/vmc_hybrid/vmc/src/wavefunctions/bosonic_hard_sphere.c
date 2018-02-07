#include "wavefunction.h"

parameters_t *get_variational_parameters(void)
{
    return (parameters_t *) 0;
}

void free_parameters_struct(parameters_t *parameters)
{
}

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
    return 0.0;
}
