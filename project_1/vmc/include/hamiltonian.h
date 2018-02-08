/*
 * Abstract class to be overwritten
 * 
 */

#pragma once

class Hamiltonian {
    public:
        Hamiltonian(class System* system);
        virtual double computeLocalEnergy(class particles* walkers) = 0;

    private:
        class System* m_system = nullptr;
};