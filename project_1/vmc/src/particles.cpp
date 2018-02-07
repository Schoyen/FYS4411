/*
    Particle/walker class
*/

#include "../include/particles.h"

// Default constructor that makes one walker in 1D
Particles::Particles() 
    : Particles(1, 1) {
} 

// Constructor for creating a certain number of walkers in d dimensions
Particles::Particles(int dimensions, int number_of_walkers) {
    m_walkers = MatrixDN::Random(dimensions, number_of_walkers);
}

// Constructor that creates n walkers in d dimensions with a given distribution spread.
Particles::Particles(int dimensions, int number_of_walkers, double distribution_spread) {
    m_walkers = MatrixDN::Random(dimensions, number_of_walkers);
    m_distribution_spread = distribution_spread;
    m_walkers = m_walkers * m_distribution_spread;
}

// Method to change the distribution spread of walkers
void Particles::setDistributionSpread(double new_spread) {
    // De-spreading with old spread, applying new spread
    m_walkers = m_walkers / m_distribution_spread;
    m_walkers = m_walkers * new_spread;
    // Storing new spread
    m_distribution_spread = new_spread;
}

// Destructor
Particles::~Particles() {
    //delete m_distribution_spread;
}