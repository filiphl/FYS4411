#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"







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
    for (int i=0; i< m_system->getNumberOfParticles(); i++){
        double ri2 = 0;
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            double* particlei = particles[i]->getPosition()[i];
            ri2 += particlei[j]*particlei[j];
        }
        argument += -m_alpha * ri2;
    }
    return exp(psi);
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
    if (m_system->analytical){
        double ddr = 0;
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                double rj = particles[i][j]*particles[i][j];
                ddr += 2*m_alpha*(2*m_alpha*rj*rj);
            }
        }
    }

    else{
        double ddr = 0;
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                double psi      =   evaluate( particles );
                particles[i]->adjustPosition( m_system -> m_stepLength, j );
                double psiPlus  =   evaluate( particles );
                particles[i]->adjustPosition( -2* m_system -> m_stepLength, j );
                double psiMinus =   evaluate( particles );
                particles[i]->adjustPosition( m_system -> m_stepLength, j );
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr/((m_system->m_stepLength)*(m_system->m_stepLength));
    }
    return ddr;
}
