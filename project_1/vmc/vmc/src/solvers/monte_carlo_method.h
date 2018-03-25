#pragma once

#include <random>

#include "wavefunction.h"

class MonteCarloMethod
{
    private:
        std::mt19937_64 m_engine;

    public:
        MonteCarloMethod();
        MonteCarloMethod(int seed);

        virtual bool step(Wavefunction *wavefunction, double step_length) = 0;

        double next_uniform()
        {
            std::uniform_real_distribution<double> dist(0, 1);
            return dist(m_engine);
        }

        double next_uniform(const double &min, const double &max)
        {
            std::uniform_real_distribution<double> dist(min, max);
            return dist(m_engine);
        }

        int next_int(const int &min, const int &max)
        {
            std::uniform_int_distribution<int> dist(min, max);
            return dist(m_engine);
        }

        double next_gaussian(const double &mean, const double &std)
        {
            std::normal_distribution<double> dist(mean, std);
            return dist(m_engine);
        }
};
