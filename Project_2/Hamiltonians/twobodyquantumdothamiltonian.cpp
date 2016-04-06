#include "twobodyquantumdothamiltonian.h"

TwoBodyQuantumDotHamiltonian::TwoBodyQuantumDotHamiltonian(System* system) :
    Hamiltonian(system)
{

}

double TwoBodyQuantumDotHamiltonian::computeLocalEnergy(std::vector<Particle *> particles)
{

    double kinetic = computeKineticEnergy(particles);

    double r2 = 0;
    for (int i=0; i<2; i++){    // Only two particles
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            r2 += particles[i]->getPosition()[j] * particles[i]->getPosition()[j];

        }
    }

    double r12 = 0;
    r12 = sqrt( (particles[0]->getPosition()[0] - particles[1]->getPosition()[0])*
            (particles[0]->getPosition()[0] - particles[1]->getPosition()[0])+
            (particles[0]->getPosition()[1] - particles[1]->getPosition()[1])*
            (particles[0]->getPosition()[1] - particles[1]->getPosition()[1]) );

    return kinetic + 0.5*r2 + 1/r12;

}

double TwoBodyQuantumDotHamiltonian::computeAnalyticalEnergy(std::vector<Particle *> particles)
{

}

