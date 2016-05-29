#ifndef MANYBODYQUANTUMDOTWAVEFUNCTION_H
#define MANYBODYQUANTUMDOTWAVEFUNCTION_H
#include "wavefunction.h"
#include <armadillo>
using namespace arma;

class ManyBodyQuantumDotWaveFunction : public WaveFunction
{

private:
    double m_oa     = 0;    // omega*alpha
    double m_wa     = 0;    // sqrt(omega*alpha)
    double m_R      = 1;
    double m_RSD    = 1;
    int    m_np     = 0;
    int    m_npHalf = 0;
    double m_derivativeStepLength = 0;
    mat    m_a;
    mat    m_slaterUp;
    mat    m_slaterDown;
    mat    m_slaterUpInverse;
    mat    m_slaterDownInverse;
    mat    m_quantumNumbers;
    double m_dUp = 0;
    double m_dDown = 0;
    double m_ddUp = 0;
    double m_ddDown = 0;
    double m_padeJastrow = 0;
    double m_dPadeJastrow = 0;

public:
    ManyBodyQuantumDotWaveFunction(class System* system, double alpha, double omega, int a, double beta);

    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    double computeGradient(std::vector<Particle *> &particles, int particle, int dimension);
    double computeRatio(std::vector<class Particle*> particles, int i, int j, double change);
    double SingleParticleWF(int nx, int ny, double x, double y);
    double ddSingleParticleWF(int i,int j);
    double slaterGrad(std::vector<class Particle*> &particles, int k, int j);
    double correlationGrad(std::vector<Particle *> &particles, int k, int j);
    double correlationLap(std::vector<class Particle*> &particles, int k);
    double hermite(int energyLevel, double position);
    double hermiteDerivative(int energyLevel, double position);
    double hermiteDoubleDerivative(int energyLevel, double position);
    double hermiteDerivativeAlpha(int energyLevel, double position);

    double f(int i, int j); // see report page 2.

    void setupSlater();
    void updateSlater(int i);
    void updateDistances(int i);
    void printParameters();
};

#endif // MANYBODYQUANTUMDOTWAVEFUNCTION_H
