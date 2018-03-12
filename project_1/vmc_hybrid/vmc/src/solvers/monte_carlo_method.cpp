#include <random>

#include "monte_carlo_method.h"



MonteCarloMethod::MonteCarloMethod()
{
    std::random_device rd;
    m_engine = std::mt19937_64(rd());
}

MonteCarloMethod::MonteCarloMethod(int seed)
{
    m_engine = std::mt19937_64(seed);
}
