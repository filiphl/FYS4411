#ifndef INTERACTINSIMPLEGAUSSIAN_H
#define INTERACTINSIMPLEGAUSSIAN_H
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

class InteractinSimpleGaussian : public WaveFunction
{
private:
    double m_alpha                  = 0;
    double m_beta                   = 0;
    double m_gamma                  = 0;
    double m_numberOfParameters     = 0;
    double m_derivativeStepLength   = 0.0001;
    double m_a                      = 0.0043;
    double dr2                      = 0;
    double absdr                    = 0;
public:
    InteractinSimpleGaussian     (System* system, double alpha, double beta, double gamma);
    double evaluate              (std::vector<class Particle*> particles);
    void computeDoubleDerivative (std::vector<Particle *> particles);


    // Get / Set
    double alpha() const;
    void setalpha(double alpha);
    double beta() const;
    void setBeta(double beta);
    double gamma() const;
    void setGamma(double gamma);
};

#endif // INTERACTINSIMPLEGAUSSIAN_H
