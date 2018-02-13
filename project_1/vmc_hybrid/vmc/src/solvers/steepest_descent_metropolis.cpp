#include "steepest_descent_metropolis.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "math_macros.h"

bool SteepestDescent::step(Wavefunction *wavefunction, double gamma) 
{

    double* parameters = wavefunction->get_parameters();

    double alpha = parameters[0];

    //wavefunction->compute_alpha_derivative();

    return false;
}

double SteepestDescent::run(
        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
        double step_length, unsigned int num_samples)
{

    /*
        In general: find expectation value of energy and search for minimum
        How to do this w/o obtaining (num) derivatives wrt. variational params? 
        Seek to minimise find a s.t. Del ev(E_L(a)) = 0. Decreases fastest in
        the direction of - Del F(a). Iterative scheme:
            a_(k+1) = a_k - gamma * Del F(a_k)
        for small gamma, F(a_k) > F(a_(k+1)), should converge.

        If |F(a_k)| >= |F(a_(k+1))| we divide gamma by two.. proceed. Crude.

        Need Initial guess for a.
    */

    return 0.0;
}
