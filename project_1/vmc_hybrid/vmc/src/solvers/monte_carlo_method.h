#pragma once

#include <random>

class MonteCarloMethod
{
    protected:
        std::mt19937_64 m_engine;
        std::uniform_int_distribution<> m_random_particle;
        std::uniform_real_distribution<double> m_random_step;

    public:
        MonteCarloMethod(unsigned int num_particles);
};
