#pragma once

class System {

    public:
        void runMetropolis (int numberOfSteps);
        
        void setWalkers(class Particles* walkers);

    private:
        int m_numberOfMetropolisSteps = 0;
        class Particles* m_walkers = nullptr;
        // int m_numberOfWalkers = 0;
        // int m_numberOfDimensions = 0;
};