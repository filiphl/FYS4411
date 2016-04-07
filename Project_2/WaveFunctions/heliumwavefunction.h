#ifndef HELIUMWAVEFUNCTION_H
#define HELIUMWAVEFUNCTION_H
#include "wavefunction.h"
#include "../system.h"

class HeliumWaveFunction : public WaveFunction
{
public:
    HeliumWaveFunction(System *system, double alpha);
    double evaluate(std::vector<Particle *> particles);
    double computeLaplacian(std::vector<Particle *> particles);
    double computeGradient(Particle *particle, int dimension);
};

#endif // HELIUMWAVEFUNCTION_H
