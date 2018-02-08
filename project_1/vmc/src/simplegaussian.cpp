/*
 * Gaussian wavefunction subclass
 * 
 */

#include "../include/simplegaussian.h"
#include "../include/wavefunction.h"
#include "../include/system.h"
#include "../include/particles.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) : WaveFunction(system) {

}

double SimpleGaussian::evaluate(Particles* particles) {


    return 0;
}

double SimpleGaussian::computeDoubleDerivative(Particles* particles) {


    return 0;
}