#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>
using namespace std;





SimpleGaussian::SimpleGaussian(System* system, double alpha) :
    WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    setAlpha(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getOldPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    double argument = 0;
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();
    for (int i=0; i< numberOfParticles; i++){
        //double ri2 = 0;
        for (int j=0; j<numberOfDimensions; j++){
            argument -= particles[i]->getOldPosition()[j] * particles[i]->getOldPosition()[j];
        }
        //argument -= ri2;
    }
    psiAlpha = argument;
    return exp(m_alpha * argument);
}



double SimpleGaussian::computeLaplacian(std::vector<class Particle*> particles) {


    double ddr = 0;

    if (m_system->getAnalyticalLaplacian()){
        cout << "ALPHA: "<<m_alpha<<endl;
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            double r2 = 0;
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                double rj = particles[i]->getOldPosition()[j];
                r2 += rj*rj;
            }
            ddr += 2*m_alpha*(2*m_alpha*r2 - m_system->getNumberOfDimensions());
        }

    }

    else {
        m_derivativeStepLength = 1e-5;
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                double psi      =   evaluate( particles );
                particles[i]->adjustOldPosition( m_derivativeStepLength, j );      // +
                double psiPlus  =   evaluate( particles );
                particles[i]->adjustOldPosition( -2 * m_derivativeStepLength, j ); // -
                double psiMinus =   evaluate( particles );
                particles[i]->adjustOldPosition( m_derivativeStepLength, j );      // reset
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr / (m_derivativeStepLength * m_derivativeStepLength);
    }
    return ddr;
}



double SimpleGaussian::computeGradient(std::vector<Particle *> particles, int particle, int dimension)
{
    return -2*m_parameters[0]*particles[particle]->getOldPosition()[dimension];
}



double SimpleGaussian::computeRatio(std::vector<class Particle*> particles, int i, int j, double change)
{
    double oldPsi = evaluate(particles);
    particles[i]->adjustOldPosition(change, j);
    double newPsi = evaluate(particles);
    return newPsi*newPsi/(oldPsi*oldPsi);
}

void SimpleGaussian::updateSlater(int i)
{

}

std::vector<double> SimpleGaussian::getParameters()
{

}


