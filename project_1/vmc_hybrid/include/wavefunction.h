#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

typedef struct particle {
    double *position;
} particle_t;

typedef struct particles {
    unsigned int num_particles;
    unsigned int dimensionality;
    particle_t *particles;
} particles_t;

typedef struct parameters {
    unsigned int num_parameters;
    double *parameters;
} parameters_t;

typedef struct wavefunction {
    parameters_t *parameters;
    particles_t *particles;
} wavefunction_t;



void allocate_variational_parameters(wavefunction_t *wavefunction);
void free_variational_parameters(wavefunction_t *wavefunction);

void allocate_particles(wavefunction_t *wavefunction);
void free_particles(wavefunction_t *wavefunction);

double ratio(
        parameters_t *parameters, double new_position[DIMENSIONALITY],
        double old_position[DIMENSIONALITY]);

double wavefunction(
        parameters_t *parameters, double position[DIMENSIONALITY]);

#endif
