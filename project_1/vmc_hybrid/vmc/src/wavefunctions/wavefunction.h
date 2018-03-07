#pragma once

#include <valarray>
#include <cmath>

#include "math_macros.h"

class Wavefunction
{
    protected:
        unsigned int m_num_particles;
        unsigned int m_num_dimensions;
        unsigned int m_num_parameters;

        double m_hard_sphere_radius = 0.0;
        double m_mass;
        double m_omega;

        double *m_parameters;

        /* We use single pointer in order to allow for allocation from NumPy.
         * This array will thus be indexed as:
         *
         *      m_particles[j + i*m_num_dimensions] = m_particles[i][j],
         *
         * where i is the particle index and j is the position index.
         *
         *      Nice to know.
         *          -Seb.
         */
        double *m_particles;

    public:
        Wavefunction(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double mass,
                double omega,
                double *parameters,
                double *particles);


        double compute_position_squared_sum();

        unsigned int get_num_particles()
        {
            return m_num_particles;
        }

        unsigned int get_num_dimensions()
        {
            return m_num_dimensions;
        }

        double *get_parameters()
        {
            return m_parameters;
        }

        unsigned int get_num_parameters()
        {
            return m_num_parameters;
        }

        double get_hard_sphere_radius()
        {
            return m_hard_sphere_radius;
        }

        double get_mass()
        {
            return m_mass;
        }

        double get_frequency()
        {
            return m_omega;
        }

        void reset_particle_position(
                double position[], unsigned int particle_index)
        {
            unsigned int i;

            for (i = 0; i < m_num_dimensions; i++) {
                m_particles[i + particle_index*m_num_dimensions] = position[i];
            }
        }

        void copy_particle_position(
                double position[], unsigned int p_k)
        {
            unsigned int i;

            for (i = 0; i < m_num_dimensions; i++) {
                position[i] = m_particles[i + p_k*m_num_dimensions];
            }
        }

        double get_distance_between_particles(
                unsigned int p_i, unsigned int p_j)
        {
            double distance;
            unsigned int i;

            distance = 0.0;

            for (i = 0; i < m_num_dimensions; i++) {
                distance += SQUARE(
                        (m_particles[i + p_i*m_num_dimensions]
                         - m_particles[i + p_j*m_num_dimensions]));
            }

            return sqrt(distance);
        }

        void move_particle(
                double step, unsigned int particle_index,
                unsigned int coordinate);

        virtual double evaluate() = 0;
        virtual double compute_laplacian() = 0;
        virtual double compute_drift_force_component(double coordinate) = 0;
        virtual std::valarray<double> compute_laplacian_variational_gradient()
            = 0;
};
