#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include <iostream>
Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeKineticEnergy(std::vector<Particle*> particles)
{
    double psi = m_system->getWaveFunction()->evaluate(particles);
    double ddr = m_system->getWaveFunction()->computeLaplacian(particles);
    double kineticEnergy = -0.5 * ddr;

    if (!(m_system->getAnalyticalLaplacian())){
         kineticEnergy /= psi;
    }
    return kineticEnergy;
}

