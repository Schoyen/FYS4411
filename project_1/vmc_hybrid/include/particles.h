#ifndef PARTICLES_H
#define PARTICLES_H

#define DIMENSIONALITY 3

typedef struct particle {
    double position[DIMENSIONALITY];
} particle_t;

typedef struct particles {
    particle_t *particles;
} particles_t;


#endif
