#ifndef HELIUMHAMILTONIAN_H
#define HELIUMHAMILTONIAN_H
#include "hamiltonian.h"
#include "../system.h"

class HeliumHamiltonian : public Hamiltonian
{
public:
    HeliumHamiltonian(System *system);
    double computeLocalEnergy(std::vector<Particle*> particles);
    double computeAnalyticalEnergy(std::vector<Particle*> particles);
    double computeGradient(Particle* particle, int dimension);
};

#endif // HELIUMHAMILTONIAN_H
