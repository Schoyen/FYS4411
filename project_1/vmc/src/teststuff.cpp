#include <iostream>
#include "../include/particles.h"

int main() {

    // Create Particle object
    Particles* some_walkers = new Particles(3, 10, 10);

    std::cout << some_walkers->get_walkers() << std::endl;

    return 0;
}