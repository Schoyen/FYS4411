/*
 * Gaussian wavefunction subclass
 * 
 */

#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
    public:
        SimpleGaussian(class System* system, double alpha);
        double evaluate(class Particles* Particles);
        double computeDoubleDerivative(class Particles* particles);
};