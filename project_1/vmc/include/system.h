#pragma once

class System {

    public:
        void runMetropolis (int numberOfSteps);
        
        void setWalkers(class Particles* walkers);
        void setHamiltonian(class Hamiltonian* hamiltonian);
        void setWaveFunction(class WaveFunction* wavefunction);

    private:
        int m_numberOfMetropolisSteps = 0;
        class Particles* m_walkers = nullptr;
        class Hamiltonian* m_hamiltonian = nullptr;
        class WaveFunction* m_wavefunction = nullptr;
        // int m_numberOfWalkers = 0;
        // int m_numberOfDimensions = 0;
};