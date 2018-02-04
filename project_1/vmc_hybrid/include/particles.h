#ifndef PARTICLES_H
#define PARTICLES_H

typedef struct particle {
    double *position;
} particle_t;

typedef struct particles {
    unsigned int num_particles;
    unsigned int dimensionality;
    particle_t *particles;
} particles_t;

#endif
