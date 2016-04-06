#include "heliumhamiltonian.h"

HeliumHamiltonian::HeliumHamiltonian(System* system) :
    Hamiltonian(system)
{
    if (system->getNumberOfParticles() != 2){
        cout << "The helium atom contains 2 electrons (particles)." << endl;
        cout << "Number of particles given: "<< system->getNumberOfParticles() << endl;
        exit(0);
    }
}

double HeliumHamiltonian::computeLocalEnergy(std::vector<Particle *> particles)
{
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            r1 += particles[0]->getPosition()[j]*particles[0]->getPosition()[j];
            r2 += particles[1]->getPosition()[j]*particles[1]->getPosition()[j];
            r12 += (particles[0]->getPosition()[j] - particles[1]->getPosition()[j])*
                   (particles[0]->getPosition()[j] - particles[1]->getPosition()[j]);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    r12 = sqrt(r12);

    return computeKineticEnergy(particles) - 2*(1/r1 + 1/r2) + 1/r12;
}

double HeliumHamiltonian::computeAnalyticalEnergy(std::vector<Particle *> particles)
{

}

