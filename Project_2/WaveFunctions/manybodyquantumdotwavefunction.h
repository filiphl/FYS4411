#ifndef MANYBODYQUANTUMDOTWAVEFUNCTION_H
#define MANYBODYQUANTUMDOTWAVEFUNCTION_H
#include "wavefunction.h"

class ManyBodyQuantumDotWaveFunction : public WaveFunction
{

private:


public:
    ManyBodyQuantumDotWaveFunction(class System* system);
    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    double computeGradient(std::vector<class Particle*> particles, int particle, int dimension);
    double computeSingleParticleWF(int nx, int ny, double x, double y);
    double hermite(int energyLevel, double position);
};

#endif // MANYBODYQUANTUMDOTWAVEFUNCTION_H
