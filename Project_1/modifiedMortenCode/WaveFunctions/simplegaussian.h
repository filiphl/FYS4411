#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
private:
    char* name = "SimpleGaussian";
    double m_alpha;
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);


    double getAlpha() const{ return m_alpha; }
    void setAlpha(double alpha){ m_alpha = alpha; }

    char *getName() const{ return name; }
};
