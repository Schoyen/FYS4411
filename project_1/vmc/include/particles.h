#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixDN; 

class Particles {
    private:
        MatrixDN m_particles;
        double m_distributionSpread = 1.0;
        int m_numberOfParticles = 0 ;
        int m_numberOfDimensions = 1;

    public:
        Particles();
        Particles(int dimensions, int m_numberOfParticles);
        Particles(int dimensions, int numberOfParticles, double distributionSpread);
        ~Particles();

        void setDistributionSpread(double newSpread);

        double rSquaredOfParticleN(int particleN);
        
        // Getters
        MatrixDN getPositions() {return m_particles;};
        double getDistributionSpread() {return m_distributionSpread;};
        int getNumberOfWalkers() {return m_numberOfParticles;};
        int getNumberOfDimensions() {return m_numberOfDimensions;}
};