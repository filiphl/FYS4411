#include "heliumwavefunction.h"

HeliumWaveFunction::HeliumWaveFunction(System* system, double alpha) :
    WaveFunction(system)
{
    if (system->getNumberOfParticles() != 2){
        cout << "The helium atom contains 2 electrons (particles)." << endl;
        cout << "Number of particles given: "<< system->getNumberOfParticles() << endl;
        exit(0);
    }
    setAlpha(alpha);
    name = "Helium wavefunction";
}

double HeliumWaveFunction::evaluate(std::vector<Particle *> particles)
{

    // exp(-a(r1+r2))
    double argument = 0;
    double r1 = 0;
    double r2 = 0;
    for (int j=0; j<m_system->getNumberOfDimensions(); j++){
        r1 += particles[0]->getOldPosition()[j]*particles[0]->getOldPosition()[j];
        r2 += particles[1]->getOldPosition()[j]*particles[1]->getOldPosition()[j];
    }
    argument += sqrt(r1) + sqrt(r2);

    return exp(-m_alpha*argument);
}

double HeliumWaveFunction::computeLaplacian(std::vector<Particle *> particles){
    if (m_system->getAnalyticalLaplacian()){
        double r1 = 0;
        double r2 = 0;
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            r1 += particles[0]->getOldPosition()[j]*particles[0]->getOldPosition()[j];
            r2 += particles[1]->getOldPosition()[j]*particles[1]->getOldPosition()[j];
        }
        r1 = sqrt(r1);
        r2 = sqrt(r2);

        return 2*m_alpha2 - 2*m_alpha/r1 - 2*m_alpha/r2;
    }
    else {
        double m_derivativeStepLength = 0.0001;         // This shouldn't be here.
        double ddr = 0;
        double psi      =   evaluate( particles );
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                particles[i]->adjustOldPosition( m_derivativeStepLength, j );      // +
                double psiPlus  =   evaluate( particles );
                particles[i]->adjustOldPosition( -2 * m_derivativeStepLength, j ); // -
                double psiMinus =   evaluate( particles );
                particles[i]->adjustOldPosition( m_derivativeStepLength, j );      // reset
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr / (m_derivativeStepLength * m_derivativeStepLength);
        return ddr;
    }
}

double HeliumWaveFunction::computeGradient(Particle *particle, int dimension)
{
    cout << "Importance sampling has not been implemented for this system."<<endl;
    exit(0);
}

double HeliumWaveFunction::computeRatio(std::vector<Particle *> particles, int i, int j, double change)
{
    double oldPsi = evaluate(particles);
    particles[i]->adjustOldPosition(change, j);
    double newPsi = evaluate(particles);
    return newPsi*newPsi/(oldPsi*oldPsi);
}

std::vector<double> HeliumWaveFunction::getParameters()
{

}

void HeliumWaveFunction::updateSlater(int i)
{

}

void HeliumWaveFunction::printParameters()
{

}


