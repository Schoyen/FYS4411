/*
    Particle/walker class
*/

#include "../include/particles.h"

// Default constructor that makes one walker in 1D
Particles::Particles() 
    : Particles(1, 1) {
} 

// Constructor for creating a certain number of walkers in d dimensions
Particles::Particles(int dimensions, int numberOfWalkers) {
    m_walkers = MatrixDN::Random(dimensions, numberOfWalkers);
    m_numberOfWalkers = numberOfWalkers;
    m_numberOfDimensions = dimensions;
}

// Constructor that creates n walkers in d dimensions with a given distribution spread.
Particles::Particles(int dimensions, int numberOfWalkers, double distributionSpread) {
    m_walkers = MatrixDN::Random(dimensions, numberOfWalkers);
    m_distributionSpread = distributionSpread;
    m_walkers = m_walkers * m_distributionSpread;
    m_numberOfWalkers = numberOfWalkers;
    m_numberOfDimensions = dimensions;
}

// Method to change the distribution spread of walkers
void Particles::setDistributionSpread(double newSpread) {
    // De-spreading with old spread, applying new spread
    m_walkers = m_walkers / m_distributionSpread;
    m_walkers = m_walkers * newSpread;
    // Storing new spread
    m_distributionSpread = newSpread;
}

// Destructor
Particles::~Particles() {
    //delete m_distribution_spread;
}