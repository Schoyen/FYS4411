#include "metropolis_algorithm.h"
#include "steepest_descent_metropolis.h"
#include "hamiltonian.h"
#include "math_macros.h"
#include <iostream>

double SteepestDescentMetropolis::steepest_descent(
        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
        double gamma, unsigned int num_samples)
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

    const int MAX_ITER = 100;

    double* parameters = wavefunction->get_parameters();

    double energy, energy_prev; 
    double gradient, gradient_prev;
    double alpha, alpha_prev;

    // It might be necessary to have gamma AND step_length
    // Metropolis stesps
    energy = run(wavefunction, hamiltonian, gamma, num_samples);

    energy_prev = energy;

    // Compute derivative of energy wrt alpha (computed analytically)
    gradient_prev = 2*(m_psiAlphaDerivative1/num_samples - (m_psiAlphaDerivative2/num_samples)*energy);

    std::cout << gradient_prev << std::endl;

    // Starting alpha (alpha-naught) 
    alpha_prev = parameters[0];

    unsigned int iteration = 0;

    while (iteration <= MAX_ITER) 
    {
        iteration++;

        // The main attraction
        alpha = alpha - gamma * gradient_prev;

        std::cout << alpha << std::endl;
    }

    // For each step I need to multiply LOCAL ENERGY with dPsi/dalpha. Summarize.
    // Must happen in metropolis loop??

    return 0.0;
}

double MetropolisAlgorithm::run(
        Wavefunction *wavefunction, Hamiltonian *hamiltonian,
        double step_length, unsigned int num_samples)
{
    double energy, local_energy;
    unsigned int i, num_accepted_states;

    /* Set initial number of accepted states */
    num_accepted_states = 0;

    /* Set initial energy */
    energy = 0;

    /* Compute initial local energy */
    local_energy = hamiltonian->compute_local_energy(wavefunction);

    /* Used by Steepest decent */
    m_psiAlphaDerivative1 = 0;
    m_psiAlphaDerivative2 = 0;

    /* Perform num_samples metropolis steps */
    for (i = 0; i < num_samples; i++) {

        /* Do a step and check if it got accepted */
        if (step(wavefunction, step_length)) {
            /* Compute new local energy */
            local_energy = hamiltonian->compute_local_energy(wavefunction);
            num_accepted_states++;
        }

        /* Add local energy */
        energy += local_energy;

        /* Adding to parameter derivatives (Steepest descent) 
           Consider making boolean so that these computations are
           only done when necessary.. */
        m_psiAlphaDerivative1 += -wavefunction->compute_position_squared_sum()*local_energy;
        m_psiAlphaDerivative2 += -wavefunction->compute_position_squared_sum();
    }

    /* Return the total energy (without normalization) */
    return energy;
}
