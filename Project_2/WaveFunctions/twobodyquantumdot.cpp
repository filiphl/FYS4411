#include "twobodyquantumdot.h"

TwoBodyQuantumDot::TwoBodyQuantumDot(System* system, double alpha, double beta, double C) :
    WaveFunction(system)
{
    if (system->getNumberOfParticles() != 2){
        cout << "The two body quantum dot contains 2 electrons (particles)." << endl;
        cout << "Number of particles given: "<< system->getNumberOfParticles() << endl;
        exit(0);
    }
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    setAlpha(alpha);
    setBeta(beta);
    setAlpha(alpha);
    m_C = C;
}

double TwoBodyQuantumDot::evaluate(std::vector<Particle *> particles)
{
    double argument1 = 0;
    for (int i=0; i<2; i++){    // Only two particles
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            argument1 += particles[i]->getPosition()[j] * particles[i]->getPosition()[j];

        }
    }
    double r12 = 0;
    r12 = sqrt( (particles[0]->getPosition()[0] - particles[1]->getPosition()[0])*
            (particles[0]->getPosition()[0] - particles[1]->getPosition()[0])+
            (particles[0]->getPosition()[1] - particles[1]->getPosition()[1])*
            (particles[0]->getPosition()[1] - particles[1]->getPosition()[1]) );

    return m_C*exp(-m_alpha*argument1*0.5)*exp(m_a*r12/(1+m_beta*r12));
}

double TwoBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles)
{
    double ddr = 0;
    if (m_system->getAnalyticalLaplacian()){
        cout << "We have not yet implemented the analytically computed Laplacian." << endl;
    }

    else {
        const double psi = evaluate( particles );

        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                particles[i]->adjustPosition( m_derivativeStepLength, j );      // +
                const double psiPlus  =   evaluate( particles );
                particles[i]->adjustPosition( -2 * m_derivativeStepLength, j ); // -
                const double psiMinus =   evaluate( particles );
                particles[i]->adjustPosition( m_derivativeStepLength, j );      // reset
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr / (m_derivativeStepLength * m_derivativeStepLength);
    }
}


double TwoBodyQuantumDot::computeDerivativeOfAlpha(){
    double sum = 0;
    cout << "YES"<<endl;
    for (Particle* particle : m_system->getParticles()){
        for (int i=0; i< m_system->getNumberOfDimensions(); i++){
            sum -= particle->getPosition()[i]*particle->getPosition()[i];

        }
    }
    return sum;
}


