#include <random>

#include "monte_carlo_method.h"

MonteCarloMethod::MonteCarloMethod(unsigned int num_particles)
{
    std::random_device rd;
    m_engine = std::mt19937(rd());

    m_random_particle = std::uniform_int_distribution<>(0, num_particles - 1);
    m_random_step = std::uniform_real_distribution<double>(0.0, 1.0);
}
