#pragma once

class Wavefunction
{
    protected:
        unsigned int m_num_particles;
        unsigned int m_dimensions;
        unsigned int m_num_parameters;

        /* This refers to the last value from the evaluation function */
        bool m_valid_last_value;
        double m_last_value;

        /* This allows us to compute the position squared once per change of
         * particle movement */
        bool m_valid_position_squared_sum;
        double m_last_position_squared_sum;

        double *m_parameters;

        /* We use single pointer in order to allow for allocation from NumPy.
         * This array will thus be indexed as:
         *
         *      m_particles[j + i*m_dimensions] = m_particles[i][j],
         *
         * where i is the particle index and j is the position index. */
        double *m_particles;

    public:
        Wavefunction(
                unsigned int num_particles,
                unsigned int dimensions,
                unsigned int num_parameters,
                double *parameters,
                double *particles);


        double compute_position_squared_sum();

        virtual double evaluate() = 0;
        virtual double compute_laplacian() = 0;
};
