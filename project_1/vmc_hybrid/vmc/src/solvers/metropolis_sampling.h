#ifndef METROPOLIS_SAMPLING_H
#define METROPOLIS_SAMPLING_H

#include <stdlib.h>
#include <limits.h>

#include "wavefunction.h"
#include "vmc_macros.h"

double perform_metropolis_step(
        wavefunction_t *wavefunction, double step_length);

double metropolis_sampling(
        wavefunction_t *wavefunction, double step_length,
        unsigned int num_samples);

#endif
