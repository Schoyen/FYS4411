#ifndef PARTICLES_H
#define PARTICLES_H

enum {
    DIMENSIONALITY = 3
};

struct particle {
    double position[DIMENSIONALITY];
};



struct particles {
    struct particle *particles;
};



#endif
