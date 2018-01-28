#ifndef PARTICLES_H
#define PARTICLES_H

enum {
    DIMENSIONALITY = 3
};

typedef struct particle {
    double position[DIMENSIONALITY];
} particle_t;



typedef struct particles {
    unsigned int num_particles;
    particle_t *particles;
} particles_t;



#endif
