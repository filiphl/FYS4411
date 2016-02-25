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
     * the particles are accessible through the particle[i].getPosition()
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
            argument -= particles[i]->getPosition()[j] * particles[i]->getPosition()[j];
        }
        //argument -= ri2;
    }
    return exp(m_alpha * argument);
}



double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    double ddr = 0;

    if (m_system->getAnalyticalDoublederivative()){

        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            double r2 = 0;
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                double rj = particles[i]->getPosition()[j];
                r2 += rj*rj;
            }
            ddr += 2*m_alpha*(2*m_alpha*r2 - m_system->getNumberOfDimensions());
        }

    }

    else {
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                double psi      =   evaluate( particles );
                particles[i]->adjustPosition( m_derivativeStepLength, j );      // +
                double psiPlus  =   evaluate( particles );
                particles[i]->adjustPosition( -2 * m_derivativeStepLength, j ); // -
                double psiMinus =   evaluate( particles );
                particles[i]->adjustPosition( m_derivativeStepLength, j );      // reset
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr / (m_derivativeStepLength * m_derivativeStepLength);
    }
    return ddr;
}
