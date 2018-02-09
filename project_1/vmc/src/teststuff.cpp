#include <iostream>
#include <vector>
#include "../include/particles.h"
#include "../include/system.h"
#include "../include/simplegaussian.h"

int main() {


    // Creating particles class instance
    int N = 5; // No of particles
    int d = 3; // No of dimensions
    Particles* particles = new Particles(d, N, 10);

    for (int  i = 0; i < N; i ++) {
        std::cout << "r^2 of particle " << i << ": " << particles->rSquaredOfParticleN(i) << std::endl;
    }

    // Creating new system instance
    System* system = new System();

    SimpleGaussian* simpleGauss = new SimpleGaussian(system, 1.0);

    std::cout << "Evauluation of simple Gaussian wavefunction: " << simpleGauss->evaluate(particles)  << std::endl;

    return 0;
}