/*
 * WaveFunction abstract superclass
 *
 * The methods evaluate and computeDoubleDerivatives are very important 
 */

#pragma once

class WaveFunction {
    public:
        WaveFunction(class System* system);

        virtual double evaluate(class Particles* particles) = 0;
        virtual double computeDoubleDerivative(class Particles* particles) = 0;

    protected:
        class System* m_system = nullptr;
};