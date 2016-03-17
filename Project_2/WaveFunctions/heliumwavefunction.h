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
};

#endif // HELIUMWAVEFUNCTION_H
