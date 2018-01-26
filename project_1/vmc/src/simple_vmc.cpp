/*
    Simple Variational Monte Carlo in 1D modelled after Thijssen.
*/

#include <iostream>
#include <math.h>
#include <armadillo>
#include <cstdlib>
#include <fstream>

// Computes sum of local energy across vector x of workers, and value of alpha.
double harmonic_osc_local_energy_sum (double alpha, arma::vec x) {

    double energy = 0;
    for (int i = 0; i < x.size(); i++) {

        energy += alpha * x(i)*x(i) *(0.5 - 2* alpha*alpha);

    }

    return energy; 

}

// Computes local energy for given x and alpha 
double harmonic_osc_local_energy (double alpha, double x) {

    return alpha * x*x *(0.5 - 2* alpha*alpha);

}

// Computes value of trial wavefunction for given paramater alpha and x
double harmonic_osc_trial_wfn (double alpha, double x) {

    return  exp(-alpha * x*x);
}

// Computes value of exact wavefunction for given x.
double harmonic_osc_exact_wfn (double x) {

    return exp(-x*x / 2);

}

// Metropolis-Hastings alogrithm function.
double metropolis (arma::vec x, int no_samples, double alpha, double step=0.05) {

    // Initiating energy change variable.
    double delta_energy = 0;

    // Seedig RNG
    srand(time(NULL));

    // Creating empty array for sampling from worker array.
    arma::vec sampling_array = arma::zeros(no_samples);

    // Filling sampling_array with possible indices of x-array (workers)
    for (int i = 0; i < no_samples; i++) {
        sampling_array(i) = (int) rand() % x.size();
    } 


    // Sampling.
    for (int i = 0; i < no_samples; i ++) {
        double random_number_1 = ((double) rand() / (RAND_MAX));
        double random_number_2 = ((double) rand() / (RAND_MAX));

        // Random x position from sampling array
        int j = sampling_array(i);

        // Changing position slightly.
        double delta_x = step * 2 * (random_number_1 - 0.5);
        double x_new = x(j) + delta_x;

        // The Metropolis-Hastings test.
        double ratio = exp(-alpha*(x_new*x_new - x(j)*x(j)));
        if (ratio >= random_number_2) {
            delta_energy = delta_energy - harmonic_osc_local_energy(alpha, x(j));
            x(j) = x_new;
            delta_energy = delta_energy + harmonic_osc_local_energy(alpha, x_new);
        }
    }

    return delta_energy;
}

int main() {

    // alpha space to scan.
    arma::vec alphas = arma::linspace(0.1, 1.6, 11); 
    int N = alphas.size();

    // Samples to take.
    int number_of_samples = 1e6;
    int number_of_walkers = 20;

    // Empty array to store energies.
    arma::vec metropolis_energies = arma::zeros(N);

    // Iterating over alpha space.
    for (int i = 0; i < N; i++) {

        // Randomising RNG seed
        arma::arma_rng::set_seed_random();
        
        // Initial walker positions.
        arma::vec x; 
        x.randu(number_of_walkers);
        x.for_each( [](arma::mat::elem_type& val) {val -= 0.5;} );

        // double initial_energy = harmonic_osc_local_energy_sum(alphas(i), x);
        double metropolis_energy = metropolis(x, number_of_samples, alphas(i));
        metropolis_energies(i) = metropolis_energy;

    }

    // Writing  energies to file.
    std::ofstream a_file;
    a_file.open ("data.dat");
    a_file << metropolis_energies;
    a_file.close();

    // Printing energies to terminal
    std::cout << metropolis_energies;

    return 0; 
}
