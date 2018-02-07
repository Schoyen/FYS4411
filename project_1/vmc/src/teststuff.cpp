#include <iostream>
#include "../include/particles.h"
#include "../include/system.h"

int main() {

    // Create Particle object
    Particles* some_walkers = new Particles(3, 2);

    // Initiate a new system
    System* system = new System();
    system.setWalkers(some_walkers);

    std::cout << some_walkers->get_walkers() << std::endl;

    return 0;
}