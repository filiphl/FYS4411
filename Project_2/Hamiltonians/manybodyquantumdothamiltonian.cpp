#include "manybodyquantumdothamiltonian.h"


ManyBodyQuantumDotHamiltonian::ManyBodyQuantumDotHamiltonian(System *system, double omega) :
    Hamiltonian(system)
{
    m_omega  = omega;
    m_omega2 = omega*omega;
}

double ManyBodyQuantumDotHamiltonian::computeLocalEnergy(std::vector<Particle *> particles)
{
    double potentialEnergy = 0;
    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        double ri2 = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            ri2 += particles[i]->getNewPosition()[d]*particles[i]->getNewPosition()[d];
        }
        potentialEnergy += 0.5*m_omega2*ri2;
        for (int j=i+1; j<m_system->getNumberOfParticles(); j++){
            potentialEnergy += 1/m_system->getWaveFunction()->m_distances(i,j);
        }
    }
    //cout << potentialEnergy << "    "<<computeKineticEnergy(particles)<<endl;
    return potentialEnergy + computeKineticEnergy(particles);
}

double ManyBodyQuantumDotHamiltonian::computeAnalyticalEnergy(std::vector<Particle *> particles)
{

}
