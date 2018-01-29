#ifndef METROPOLIS_SAMPLING_H
#define METROPOLIS_SAMPLING_H

#include <stdlib.h>
#include <limits.h>

#include "particles.h"

/* Create a macro for creating a uniform double in the interval [0, 1) */
#define RANDOM_UNIFORM_DOUBLE \
    (((((unsigned long) arc4random() << 32) | arc4random()) \
     / ((double) UINT64_MAX)))

double perform_metropolis_step(
        particles_t *particles, double step_length,
        double (*local_energy)(double position),
        double (*ratio)(double new_position, double old_position));
void metropolis_sampling(void);

#endif
