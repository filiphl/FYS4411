#ifndef TWOBODYQUANTUMDOTHAMILTONIAN_H
#define TWOBODYQUANTUMDOTHAMILTONIAN_H
#include "hamiltonian.h"
#include "../system.h"

class TwoBodyQuantumDotHamiltonian : public Hamiltonian
{
private:
    double m_omega = 0;
public:
    TwoBodyQuantumDotHamiltonian(System *system, double omega);
    double computeLocalEnergy(std::vector<class Particle*> particles) ;
    double computeAnalyticalEnergy(std::vector<class Particle*> particles);
};

#endif // TWOBODYQUANTUMDOTHAMILTONIAN_H
