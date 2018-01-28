#ifndef METROPOLIS_SAMPLING_H
#define METROPOLIS_SAMPLING_H

#include <stdlib.h>
#include <limits.h>

#define RANDOM_UNIFORM_DOUBLE \
    (((((unsigned long) arc4random() << 32) | arc4random()) \
     / ((double) UINT64_MAX)))

double perform_metropolis_step(double step_length);
void metropolis_sampling();

#endif
