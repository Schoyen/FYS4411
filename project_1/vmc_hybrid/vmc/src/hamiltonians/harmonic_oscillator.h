#ifndef HARMONIC_OSCILLATOR_H
#define HARMONIC_OSCILLATOR_H

#include "hamiltonian.h"

typedef struct hamiltonian {
    double hbar;
    double mass;
    double omega;
} hamiltonian_t;

#endif
