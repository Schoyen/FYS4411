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
        double **m_particles;

    public:
        Wavefunction(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double mass,
                double omega,
                double *parameters,
                double *particles);

        virtual ~Wavefunction();


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

        double get_radial_position(unsigned int p_i)
        {
            unsigned int i;
            double radius = 0;

            for (i = 0; i < m_num_dimensions; i++) {
                radius += SQUARE(m_particles[p_i][i]);
            }

            return sqrt(radius);
        }

        void reset_particle_position(double *position, unsigned int p_i)
        {
            unsigned int i;

            for (i = 0; i < m_num_dimensions; i++) {
                m_particles[p_i][i] = position[i];
            }
        }

        void copy_particle_position(double *position, unsigned int p_k)
        {
            unsigned int i;

            for (i = 0; i < m_num_dimensions; i++) {
                position[i] = m_particles[p_k][i];
            }
        }

        double get_distance_between_particles(
                unsigned int p_i, unsigned int p_j)
        {
            double distance;
            unsigned int i;

            distance = 0.0;

            for (i = 0; i < m_num_dimensions; i++) {
                distance += SQUARE(m_particles[p_i][i] - m_particles[p_j][i]);
            }

            return sqrt(distance);
        }

        void move_particle(double step, unsigned int p_i, unsigned int i)
        {
            m_particles[p_i][i] += step;
        }

        double compute_drift_force_component(unsigned int p_i, unsigned int i)
        {
            return 2*compute_gradient_component(p_i, i);
        }

        void compute_drift_force(double *drift_force, unsigned int p_i)
        {
            unsigned int i;

            for (i = 0; i < m_num_dimensions; i++) {
                drift_force[i] = compute_drift_force_component(p_i, i);
            }
        }

        virtual double evaluate() = 0;
        virtual double compute_laplacian() = 0;
        virtual double compute_gradient_component(
                unsigned int p_i, unsigned int i) = 0;
        virtual std::valarray<double> compute_variational_gradient() = 0;
};
