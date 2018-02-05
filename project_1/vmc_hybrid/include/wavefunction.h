#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#define SQUARE(x) (x*x)

typedef struct wavefunction {
    unsigned int num_particles;
    unsigned int dimensionality;
    unsigned int num_parameters;

    double last_value;

    double *parameters;
    double **particles;
} wavefunction_t;



void allocate_variational_parameters(wavefunction_t *wavefunction);
void free_variational_parameters(wavefunction_t *wavefunction);

void allocate_particles(wavefunction_t *wavefunction);
void free_particles(wavefunction_t *wavefunction);

double local_energy(wavefunction_t *wavefunction);

double ratio(wavefunction_t *wavefunction);

double evaluate_wavefunction(wavefunction_t *wavefunction);

#endif
