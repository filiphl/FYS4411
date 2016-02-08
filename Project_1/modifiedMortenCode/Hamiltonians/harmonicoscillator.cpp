#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <iomanip>

using std::cout;
using std::endl;
using namespace std;


HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double potentialEnergy = 0;
    double kineticEnergy   = 0;

    // Compute potential energy
    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        std::vector<double> position = particles.at(i)->getPosition();  // vector
        double r2 = 0;
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            r2 += position[j]*position[j];
        }
        potentialEnergy += 0.5*m_omega*m_omega*r2;
    }

    kineticEnergy = computeKineticEnergy(particles);
    /*
    cout << "kinetic: " << setw(10) << setprecision(3) << kineticEnergy;
    cout <<"   potential: " <<  setw(10) << setprecision(3) << potentialEnergy;
    cout <<"   total: " << setw(10) << setprecision(3) << kineticEnergy + potentialEnergy << endl;
    */
    return kineticEnergy + potentialEnergy;
}

