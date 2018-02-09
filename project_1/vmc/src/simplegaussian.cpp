/*
 * Gaussian wavefunction subclass
 * 
 */

#include "../include/simplegaussian.h"
#include "../include/wavefunction.h"
#include "../include/system.h"
#include "../include/particles.h"
#include <cmath>
#include <iostream>

SimpleGaussian::SimpleGaussian(System* system, double alpha) : WaveFunction(system) {
    m_numberOfParameters = 1;
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(Particles* particles) {

    // Not sure I need this due to the nature of getting r**2
    // Assuming we have only one particle to begin with?
    // int numberOfDimensions = particles->getNumberOfDimensions();
    
     // The deal: psi = exp(-alpha * r * r)
    int numberOfParticles = particles->getNumberOfParticles();
    std::cout << m_parameters[1] << std::endl;
    std::cout << "Size : " << m_parameters.size() << std::endl;

    double psi; 
    for (int i = 0; i < numberOfParticles; i++) {
        double rSquared = particles->rSquaredOfParticleN(i);
        psi += std::exp(-m_parameters[0]*rSquared);
        std::cout << psi << std::endl;
    }

    return psi;
}

double SimpleGaussian::computeDoubleDerivative(Particles* particles) {


    return 0;
}