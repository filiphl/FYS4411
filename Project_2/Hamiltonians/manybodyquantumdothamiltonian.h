#ifndef MANYBODYQUANTUMDOTHAMILTONIAN_H
#define MANYBODYQUANTUMDOTHAMILTONIAN_H
#include "hamiltonian.h"

class ManyBodyQuantumDotHamiltonian : public Hamiltonian
{
private:
    double m_omega  = 0;
    double m_omega2 = 0;    // These must be set from main.
public:
    ManyBodyQuantumDotHamiltonian(class System* system, double omega);
    double computeLocalEnergy(std::vector<class Particle*> particles);
    double computeAnalyticalEnergy(std::vector<class Particle*> particles);
};

#endif // MANYBODYQUANTUMDOTHAMILTONIAN_H
