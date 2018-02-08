/*
 * Gaussian wavefunction subclass
 * 
 */

#include "../include/simplegaussian.h"
#include "../include/wavefunction.h"
#include "../include/system.h"
#include "../include/particles.h"
#include <cmath>

SimpleGaussian::SimpleGaussian(System* system, double alpha) : WaveFunction(system) {

}

double SimpleGaussian::evaluate(Particles* particles) {

    // Assuming we have only one particle to begin with?
    int numberOfDimensions = particles->getNumberOfDimensions;

    double r_squared = 0;

    // The deal: psi = exp(-alpha * r * r)
    return 0;
}

double SimpleGaussian::computeDoubleDerivative(Particles* particles) {


    return 0;
}