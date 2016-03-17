#ifndef INTERACTINGHARMONICOSCILLATOR_H
#define INTERACTINGHARMONICOSCILLATOR_H
#include "hamiltonian.h"

class InteractingHarmonicOscillator : public Hamiltonian
{
private:
    double m_omegaHO2 = 0;
    double m_omegaZ2  = 0;
    double m_gamma2   = 0;
    double m_a        = 0.0043;
public:
    InteractingHarmonicOscillator(System *system, double omegaHO, double omegaZ, double gamma);
    double computeLocalEnergy(std::vector<Particle*> particles);
    double computeAnalyticalEnergy(std::vector<Particle*> particles);

    double omegaHO() const;
    void setOmegaHO(double omegaHO);
    double omegaZ() const;
    void setOmegaZ(double omegaZ);
    double gamma() const;
    void setGamma(double gamma);
    double a() const;
    void setA(double a);
};

#endif // INTERACTINGHARMONICOSCILLATOR_H
