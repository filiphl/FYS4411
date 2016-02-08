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
    double kineticEnergy = 0;
    for (int i=0; i< m_system->getNumberOfParticles(); i++){
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            double ddr = m_system->getWaveFunction()->computeDoubleDerivative(particles);
            kineticEnergy += ddr;
        }
    }
    kineticEnergy *= -0.5 / psi;
    return kineticEnergy;
}

