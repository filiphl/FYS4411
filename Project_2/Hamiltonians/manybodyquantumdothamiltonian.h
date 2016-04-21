#ifndef MANYBODYQUANTUMDOTHAMILTONIAN_H
#define MANYBODYQUANTUMDOTHAMILTONIAN_H
#include "hamiltonian.h"

class ManyBodyQuantumDotHamiltonian : public Hamiltonian
{
public:
    ManyBodyQuantumDotHamiltonian(class System* system);
    double computeLocalEnergy(std::vector<class Particle*> particles);
    double computeAnalyticalEnergy(std::vector<class Particle*> particles);
};

#endif // MANYBODYQUANTUMDOTHAMILTONIAN_H
