#include "interactingharmonicoscillator.h"





double InteractingHarmonicOscillator::a() const
{
    return m_a;
}

void InteractingHarmonicOscillator::setA(double a)
{
    m_a = a;
}

InteractingHarmonicOscillator::InteractingHarmonicOscillator(System* system, double omegaHO, double omegaZ, double gamma):
    Hamiltonian(system)
{
    assert(omegaHO > 0);
    m_omegaHO2  = omegaHO*omegaHO;
    assert(omegaZ > 0);
    m_omegaZ2   = omegaZ*omegaZ;
    m_gamma2   = gamma*gamma;
}

double InteractingHarmonicOscillator::computeLocalEnergy(std::vector<Particle *> particles)
{
    double potentialEnergy = 0;
    double kineticEnergy   = 0;

    // Compute external potential energy
    double r2 = 0;
    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            if (j<2){
                r2 += particles[i]->getPosition()[j] * particles[i]->getPosition()[j];
            }
            else{
                r2 += particles[i]->getPosition()[j] * particles[i]->getPosition()[j] * m_gamma2;
            }
        }
    }

    potentialEnergy = 0.5*r2;

    // Compute internal potential energy
    double dr2 = 0;
    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        for (int j=i+1; j<m_system->getNumberOfParticles(); j++){
            dr2=0;
            for (int k=0; k<m_system->getNumberOfDimensions(); k++){
                dr2 += (particles[j]->getPosition()[k] - particles[i]->getPosition()[k]) *
                       (particles[j]->getPosition()[k] - particles[i]->getPosition()[k]) ;
            }
            double absdr = sqrt(dr2);
            if (absdr < m_a){
                potentialEnergy += 1e10;
            }
        }
    }

    kineticEnergy = computeKineticEnergy(particles);
    /*
    cout << setw(10) << " kinetic: "   << setw(10) << setprecision(3) << left << kineticEnergy;
    cout << setw(12) << " potential: " << setw(10) << setprecision(3) << left << potentialEnergy;
    cout << setw(8) <<  " total: "     << setw(10) << setprecision(3) << left << kineticEnergy + potentialEnergy << endl;
    */

    return kineticEnergy + potentialEnergy;
}

double InteractingHarmonicOscillator::computeAnalyticalEnergy(std::vector<Particle *> particles)
{
    return 0;
}








// Get / Set

double InteractingHarmonicOscillator::omegaHO() const
{
    return m_omegaHO2;
}

void InteractingHarmonicOscillator::setOmegaHO(double omegaHO)
{
    m_omegaHO2 = omegaHO;
}

double InteractingHarmonicOscillator::omegaZ() const
{
    return m_omegaZ2;
}

void InteractingHarmonicOscillator::setOmegaZ(double omegaZ)
{
    m_omegaZ2 = omegaZ;
}

double InteractingHarmonicOscillator::gamma() const
{
    return m_gamma2;
}

void InteractingHarmonicOscillator::setGamma(double gamma)
{
    m_gamma2 = gamma;
}
