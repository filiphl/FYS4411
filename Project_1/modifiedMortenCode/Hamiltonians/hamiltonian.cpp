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
    double ddr = m_system->getWaveFunction()->computeDoubleDerivative(particles);
    double kineticEnergy = -0.5 * m_system->getWaveFunction()->computeDoubleDerivative(particles);
    if (!(m_system->getAnalyticalDoublederivative())){
        kineticEnergy /= psi;
    }
    return kineticEnergy;
}

