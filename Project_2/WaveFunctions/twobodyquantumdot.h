#ifndef TWOBODYQUANTUMDOT_H
#define TWOBODYQUANTUMDOT_H
#include "wavefunction.h"
#include "../system.h"

class TwoBodyQuantumDot : public WaveFunction
{
private:
    double m_a        = 1;
    double m_C        = 0;
    double m_derivativeStepLength = 0.001;
public:
    TwoBodyQuantumDot(System *system, double alpha, double beta, double C);
    double evaluate(std::vector<Particle *> particles);
    double computeLaplacian(std::vector<Particle *> particles);
    double computeDerivativeOfAlpha();
};

#endif // TWOBODYQUANTUMDOT_H