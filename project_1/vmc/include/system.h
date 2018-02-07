#include "./particles.h"

class System {

    public:
        void runMetropolis (int numberOfSteps);
        
        void setWalkers(class Particles* walkers);

    private:
        Particles m_walkers;
        int m_numberOfWalkers = 0;
        int m_numberOfDimensions = 0;
};