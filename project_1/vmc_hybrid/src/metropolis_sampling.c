#include "metropolis_sampling.h"

double perform_metropolis_step(double step_length)
{
    double step = step_length*(2.0*RANDOM_UNIFORM_DOUBLE - 1.0);
    return step;
}
{
}
