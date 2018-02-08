/*
 * Class that defines that system, also contains solver
 *
 */

#include "../include/system.h"


void System::runMetropolis(int numberOfSteps) {
    m_numberOfMetropolisSteps = numberOfSteps;

    for (int i = 0; i < numberOfSteps; i++) {
        
        // Trying to get an energy
        double localEnergy = 2;
    }
};

void System::setWalkers(Particles* walkers) {
    m_walkers = walkers;
};

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    // Stuff
}

