/*
    Particle/walker class
*/

#include "../include/particles.h"

// Default constructor that makes two walkers in 1D
Particles::Particles() 
    : Particles(1, 2) {
} 

Particles::Particles(int dimensions, int number_of_walkers) {
    m_walkers = MatrixDN::Random(dimensions, number_of_walkers);
}

Particles::Particles(int dimensions, int number_of_walkers, double distribution_spread) {
    m_walkers = MatrixDN::Random(dimensions, number_of_walkers);
    m_distribution_spread = distribution_spread;
    m_walkers = m_walkers * m_distribution_spread;
}

void Particles::set_distribution_spread(double new_spread) {
    // De-spreading with old spread, applying new spread
    m_walkers = m_walkers / m_distribution_spread;
    m_walkers = m_walkers * new_spread;
    // Storing new spread
    m_distribution_spread = new_spread;
}

Particles::~Particles() {
    //delete m_distribution_spread;
}