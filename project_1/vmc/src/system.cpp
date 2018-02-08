/*
 * Class that defines that system, also contains solver.
 * 
 * Needed: Hamiltonian, WaveFunction, Particles/Walkers
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

void System::setParticles(Particles* particles) {
    m_particles = particles;
};

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
};

void System::setWaveFunction(WaveFunction* wavefunction) {
    m_wavefunction = wavefunction;
};

