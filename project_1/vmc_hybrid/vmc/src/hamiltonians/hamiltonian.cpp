#include <cmath>

#include "hamiltonian.h"
#include "wavefunction.h"



double Hamiltonian::compute_local_energy(Wavefunction *wavefunction)
{
    double kinetic_energy, potential_energy;

    kinetic_energy = compute_kinetic_energy(wavefunction);
    potential_energy = compute_potential_energy(wavefunction);

    return kinetic_energy + potential_energy;
}
