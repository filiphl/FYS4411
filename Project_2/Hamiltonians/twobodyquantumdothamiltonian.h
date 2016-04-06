#ifndef TWOBODYQUANTUMDOTHAMILTONIAN_H
#define TWOBODYQUANTUMDOTHAMILTONIAN_H
#include "hamiltonian.h"
#include "../system.h"

class TwoBodyQuantumDotHamiltonian : public Hamiltonian
{
private:

public:
    TwoBodyQuantumDotHamiltonian(System *system);
    double computeLocalEnergy(std::vector<class Particle*> particles) ;
    double computeAnalyticalEnergy(std::vector<class Particle*> particles);
};

#endif // TWOBODYQUANTUMDOTHAMILTONIAN_H
