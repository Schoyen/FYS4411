#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixDN; 

class Particles {
    private:
        MatrixDN m_walkers;
        double m_distribution_spread;

    public:
        Particles();
        Particles(int dimensions, int number_of_walkers);
        Particles(int dimensions, int number_of_walkers, double distribution_spread);
        ~Particles();

        void set_distribution_spread(double new_spread);
        
        // Getters
        MatrixDN get_walkers() {return m_walkers;};
        double get_distribution_spread() {return m_distribution_spread;};
};