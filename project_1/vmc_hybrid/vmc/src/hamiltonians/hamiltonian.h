#pragma once

class Hamiltonian
{
    public:
        virtual double compute_local_energy(Wavefunction *wavefunction) = 0;
};
