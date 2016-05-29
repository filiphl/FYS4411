#ifndef TWOBODYQUANTUMDOT_H
#define TWOBODYQUANTUMDOT_H
#include "wavefunction.h"
#include "../system.h"

class TwoBodyQuantumDot : public WaveFunction
{
private:
    double m_a        = 0;
    double m_C        = 0;
    double m_derivativeStepLength = 0;
    double m_omega    = 0;
public:
    TwoBodyQuantumDot(System *system, double alpha, double beta, double C, double omega, double a);
    double evaluate(std::vector<Particle *> particles);
    double computeLaplacian(std::vector<Particle *> particles);
    double computeGradient(std::vector<Particle *>& particles, int particle, int dimension);
    double computeRatio(std::vector<class Particle*> particles, int i, int j, double change);
    double computeDerivativeOfAlpha();
    std::vector<double> getParameters();
    void   updateSlater(int i);
    void   printParameters();
};

#endif // TWOBODYQUANTUMDOT_H
