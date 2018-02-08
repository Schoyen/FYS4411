#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixDN; 

class Particles {
    private:
        MatrixDN m_walkers;
        double m_distributionSpread = 1.0;
        int m_numberOfWalkers = 0 ;
        int m_numberOfDimensions = 1;

    public:
        Particles();
        Particles(int dimensions, int numberOfWalkers);
        Particles(int dimensions, int numberOfWalkers, double distributionSpread);
        ~Particles();

        void setDistributionSpread(double newSpread);
        
        // Getters
        MatrixDN get_walkers() {return m_walkers;};
        double getDistributionSpread() {return m_distributionSpread;};
        int getNumberOfWalkers() {return m_numberOfWalkers;}
        int getNumberOfDimensions() {return m_numberOfDimensions;}
};