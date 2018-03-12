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
};
