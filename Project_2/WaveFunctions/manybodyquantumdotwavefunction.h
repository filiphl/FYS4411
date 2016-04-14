#ifndef MANYBODYQUANTUMDOTWAVEFUNCTION_H
#define MANYBODYQUANTUMDOTWAVEFUNCTION_H
#include "wavefunction.h"

class ManyBodyQuantumDotWaveFunction// : public WaveFunction
{

private:


public:
    ManyBodyQuantumDotWaveFunction(class System* system);

    double evaluate(std::vector<class Particle*> particles, int energyLevel);
    double computeLaplacian(std::vector<class Particle*> particles);
    double computeGradient(std::vector<class Particle*> particles, int particle, int dimension);
    double computeRatio(std::vector<class Particle*> particles, int i, int j, double change);

};

#endif // MANYBODYQUANTUMDOTWAVEFUNCTION_H
