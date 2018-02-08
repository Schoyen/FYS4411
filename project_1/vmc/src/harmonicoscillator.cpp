#include "../include/harmonicoscillator.h"
#include "../include/hamiltonian.h"
#include "../include/system.h"
#include "../include/particles.h"

HarmonicOscillator::HarmonicOscillator(System* system, double omega) : Harmiltonian(system) {
    m_omega = omega;
}  

double HarmonicOscillator::computeLocalEnergy() {
    
}