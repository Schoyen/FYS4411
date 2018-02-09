/*
 * WaveFunction abstract superclass
 *
 * The methods evaluate and computeDoubleDerivatives are very important 
 */

#pragma once
#include <vector>

class WaveFunction {
    public:
        WaveFunction(class System* system);

        virtual double evaluate(class Particles* particles) = 0;
        virtual double computeDoubleDerivative(class Particles* particles) = 0;

    protected:
        int m_numberOfParameters = 0;
        class System* m_system = nullptr;
        std::vector<double> m_parameters = std::vector<double>();
        
};