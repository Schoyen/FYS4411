#ifndef METROPOLIS_SAMPLING_H
#define METROPOLIS_SAMPLING_H

#include <stdlib.h>
#include <limits.h>

#include "wavefunction_config.h"
#include "particles.h"
#include "wavefunction.h"

/* Create a macro for creating a uniform double in the interval [0, 1) */
#define RANDOM_UNIFORM_DOUBLE \
    (((((unsigned long) arc4random() << 32) | arc4random()) \
     / ((double) UINT64_MAX)))

double perform_metropolis_step(
        particles_t *particles,
        parameters_t *parameters,
        double step_length);

double metropolis_sampling(
        particles_t *particles,
        parameters_t *parameters,
        double step_length,
        unsigned int num_samples);

#endif
