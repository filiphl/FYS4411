#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
private:
    //double m_alpha;
    double m_derivativeStepLength = 0;
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    double computeGradient(std::vector<class Particle*> particles, int particle, int dimension);

    //double getAlpha() const{ return m_alpha; }
    //void setAlpha(double alpha){ m_alpha = alpha; }

};
