#ifndef MANYBODYQUANTUMDOTWAVEFUNCTION_H
#define MANYBODYQUANTUMDOTWAVEFUNCTION_H
#include "wavefunction.h"
#include <armadillo>
using namespace arma;

class ManyBodyQuantumDotWaveFunction : public WaveFunction
{

private:
    double m_omega  = 0;
    double m_beta   = 0;
    double m_R      = 1;
    int    m_npHalf = 0;
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
    ManyBodyQuantumDotWaveFunction(class System* system, double omega, double a, double beta);

    double evaluate(std::vector<class Particle*> particles);
    double computeLaplacian(std::vector<class Particle*> particles);
    double computeGradient(std::vector<class Particle*> particles, int particle, int dimension);
    double computeRatio(std::vector<class Particle*> particles, int i, int j, double change);
    double SingleParticleWF(int nx, int ny, double x, double y);
    double ddSingleParticleWF(int i,int j);
    double slaterGrad(std::vector<class Particle*> particles, int k, int j);
    double correlationGrad(std::vector<class Particle*> particles, int k, int j);
    double correlationLap(std::vector<class Particle*> particles, int k);
    double hermite(int energyLevel, double position);
    double hermiteDerivative(int energyLevel, double position);
    double hermiteDoubleDerivative(int energyLevel, double position);

    double f(int i, int j); // see report page 2.


    void setupSlater();
    void updateSlater(int i);
    void printParameters();
};

#endif // MANYBODYQUANTUMDOTWAVEFUNCTION_H
