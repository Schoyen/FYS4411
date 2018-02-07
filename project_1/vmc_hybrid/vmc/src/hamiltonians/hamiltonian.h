#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "wavefunction.h"

typedef struct hamiltonian hamiltonian_t;

double local_energy(wavefunction_t *wavefunction, hamiltonian_t *hamiltonian);

#endif
