#pragma once

#include "wavefunction.h"

class SimpleGaussian : public Wavefunction
{
    public:
        using Wavefunction::Wavefunction;

        double evaluate();
        double compute_laplacian();
};
