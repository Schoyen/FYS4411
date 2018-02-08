/*
    Particle/walker class
*/

#include "../include/particles.h"

// Default constructor that makes one walker in 1D
Particles::Particles() 
    : Particles(1, 1) {
} 

// Constructor for creating a certain number of walkers in d dimensions
Particles::Particles(int dimensions, int numberOfParticles) {
    m_particles = MatrixDN::Random(dimensions, numberOfParticles);
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = dimensions;
}

// Constructor that creates n walkers in d dimensions with a given distribution spread.
Particles::Particles(int dimensions, int numberOfParticles, double distributionSpread) {
    m_particles = MatrixDN::Random(dimensions, numberOfParticles);
    m_distributionSpread = distributionSpread;
    m_particles = m_particles * m_distributionSpread;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = dimensions;
}

// Method to change the distribution spread of walkers
void Particles::setDistributionSpread(double newSpread) {
    // De-spreading with old spread, applying new spread
    m_particles = m_particles / m_distributionSpread;
    m_particles = m_particles * newSpread;
    // Storing new spread
    m_distributionSpread = newSpread;
}

// Destructor
Particles::~Particles() {
    //delete m_distribution_spread;
}