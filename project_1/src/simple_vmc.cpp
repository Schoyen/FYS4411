#include <iostream>
#include <math.h>
#include <armadillo>
#include <cstdlib>
#include <fstream>

double harmonic_osc_local_energy_sum (double alpha, arma::vec x) {

    double energy = 0;
    for (int i = 1; i <= x.size(); i++) {

        std::cout << "hei" << std::endl;
        energy += alpha * x(i)*x(i) *(0.5 - 2* alpha*alpha);

    }

    std::cout << x << std::endl;

    return energy; 

}

double harmonic_osc_local_energy (double alpha, double x) {

    return alpha * x*x *(0.5 - 2* alpha*alpha);

}

double harmonic_osc_trial_wfn (double alpha, double x) {

    return  exp(-alpha * x*x / 2);

}

double harmonic_osc_exact_wfn (double x) {

    return exp(-x*x / 2);

}

double metropolis (arma::vec x, int no_samples, double alpha, double step=0.05) {

    double delta_energy = 0;

    // Seedig RNG
    srand(time(NULL));

    arma::vec sampling_array = arma::zeros(no_samples);

    // Filling sampling_array
    for (int i = 1; i < no_samples; i++) {
        sampling_array(i) = (int) rand() % x.size();
    } 

    for (int i = 1; i < no_samples; i ++) {
        std::cout << no_samples - i << std::endl; 
        double random_number_1 = ((double) rand() / (RAND_MAX));
        double random_number_2 = ((double) rand() / (RAND_MAX));

        double delta_x = step * 2 * (random_number_1 - 0.5);
        double x_new = x(i) + delta_x;

        // The test
        double ratio = exp(-alpha*(x_new*x_new - x(i)*x(i)));
        if (ratio >= random_number_2) {
            delta_energy = delta_energy - harmonic_osc_local_energy(alpha, x(i));
            x(i) = x_new;
            delta_energy = delta_energy + harmonic_osc_local_energy(alpha, x_new);
        }
    }

    std::cout << delta_energy << std::endl;

    return delta_energy;
}

int main() {

    arma::vec alphas = arma::linspace(0.1, 1.6, 11); 
    int N = alphas.size();
    int number_of_samples = 1000000;
    int number_of_walkers = 2;

    arma::vec energies = arma::zeros(N);
    arma::vec metropolis_energies = arma::zeros(N);

    for (int i = 1; i < N; i++) {
        
        std::cout << alphas(i) << std::endl;
        arma::arma_rng::set_seed_random();
        arma::vec x; 
        x.randu(2);
        x.for_each( [](arma::mat::elem_type& val) {val -= 0.5;} );

        std::cout << x.size() << std::endl;

        double initial_energy = harmonic_osc_local_energy_sum(alphas(i), x);
        double metropolis_energy = metropolis(x, number_of_samples, alphas(i));
        metropolis_energies(i) = metropolis_energy;

    }

    std::ofstream a_file;
    a_file.open ("data.dat");

    a_file << metropolis_energies;

    a_file.close();

    return 0; 
}
