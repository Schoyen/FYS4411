#ifndef SIMPLE_VMC_H
#define SIMPLE_VMC_H

#include <armadillo>

double harmonic_osc_local_energy_sum (double alpha, arma::vec x);
double harmonic_osc_local_energy (double alpha, double x);
double harmonic_osc_trial_wfn (double alpha, double x);
double harmonic_osc_exact_wfn (double x);
double metropolis (arma::vec x, int no_samples, double alpha, double step);

#endif
