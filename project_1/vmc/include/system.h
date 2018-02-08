#pragma once
#include "particles.h"

class System {

    public:
        void runMetropolis (int numberOfSteps);
        
        void setWalkers(class Particles* walkers);

    private:
        class Particles* m_walkers = nullptr;
        int m_numberOfWalkers = 0;
        int m_numberOfDimensions = 0;
};