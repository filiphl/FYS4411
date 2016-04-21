#ifndef MANYBODYQUANTUMDOTWAVEFUNCTION_H
#define MANYBODYQUANTUMDOTWAVEFUNCTION_H
#include "wavefunction.h"
#include <armadillo>
using namespace arma;

class ManyBodyQuantumDotWaveFunction : public WaveFunction
{

private:
    double m_omega = 0;
    mat    m_slaterUp;
    mat    m_slaterDown;
public:
    ManyBodyQuantumDotWaveFunction(class System* system);

    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    double computeGradient(std::vector<class Particle*> particles, int particle, int dimension);
    double computeRatio(std::vector<class Particle*> particles, int i, int j, double change);
    double computeSingleParticleWF(int nx, int ny, double x, double y);
    double hermite(int energyLevel, double position);


    void setupSlater();
};

#endif // MANYBODYQUANTUMDOTWAVEFUNCTION_H
