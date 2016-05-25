#include "twobodyquantumdothamiltonian.h"

TwoBodyQuantumDotHamiltonian::TwoBodyQuantumDotHamiltonian(System* system, double omega) :
    Hamiltonian(system)
{
    m_omega = omega;
}

double TwoBodyQuantumDotHamiltonian::computeLocalEnergy(std::vector<Particle *> particles)
{

    double kinetic = computeKineticEnergy(particles);

    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    for (int i=0; i<m_system->getNumberOfDimensions(); i++){
        r1 += particles[0]->getOldPosition()[i]*particles[0]->getOldPosition()[i];
        r2 += particles[1]->getOldPosition()[i]*particles[1]->getOldPosition()[i];
        r12 += (particles[0]->getOldPosition()[i]-particles[1]->getOldPosition()[i])*
                (particles[0]->getOldPosition()[i]-particles[1]->getOldPosition()[i]);
    }
    r12 = sqrt(r12);
    double potentialEnergy = 0.5*m_omega*m_omega*(r1+r2);// + 1/r12;

    return kinetic + potentialEnergy;
}

double TwoBodyQuantumDotHamiltonian::computeAnalyticalEnergy(std::vector<Particle *> particles)
{

}


