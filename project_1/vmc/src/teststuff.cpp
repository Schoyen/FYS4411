#include <iostream>
#include "../include/particles.h"
#include "../include/system.h"

int main() {


    // Creating particles class instance
    int N = 5; // No of particles
    int d = 3; // No of dimensions
    Particles* particles = new Particles(d, N);

    int n = 2;

    for (int  i = 0; i < N; i ++) {
        std::cout << "r^2 of particle " << i << ": " << particles->rSquaredOfParticleN(i) << std::endl;
    }
    /*
    // Testing Eigen functionality

    Particles* particles = new Particles(3,5);
    MatrixDN thePostitions = particles->getPositions();

    // Getting vector
    Eigen::Vector3d positionN = thePostitions.col(2);

    std::cout << positionN << std::endl;
    std::cout << "Dot product is: " << positionN.dot(positionN) << std::endl;

    // Create Particle object
    Particles* some_walkers = new Particles(3, 2);

    // Initiate a new system
    System* system = new System();
    system->setWalkers(some_walkers);

    int numberOfWalkers = system->get_walkers

    */
    return 0;
}