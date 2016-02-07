#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include <iostream>
Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
    m_waveFunction = system->getWaveFunction();
}

void Hamiltonian::computeKineticEnergy(std::vector<Particle*> particles)
{
    double kinetickEnergy = 0;
    for (int i=0; i< m_system->getNumberOfParticles(); i++){
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            double ddr = m_waveFunction->computeDoubleDerivative(particles);

        }
    }
}

