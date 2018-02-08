cdef extern from "constants.h":
    cdef unsigned int HBAR

cdef extern from "wavefunction.h":
    cdef cppclass Wavefunction:
        Wavefunction(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double *parameters,
                double *particles) except +

        double compute_position_squared_sum()
        unsigned int get_num_dimensions()

        void reset_particle_position(
                double position[], unsigned int particle_index)

        void copy_particle_position(
                double position[], unsigned int particle_index)

        void add_step(
                double step, unsigned int particle_index,
                unsigned int coordinate)

cdef extern from "simple_gaussian.h":
    cdef cppclass SimpleGaussian(Wavefunction):
        SimpleGaussian(
                unsigned int num_particles,
                unsigned int num_dimensions,
                unsigned int num_parameters,
                double *parameters,
                double *particles) except +

        double evaluate()
        double compute_laplacian()
