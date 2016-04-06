#include "twobodyquantumdot.h"

TwoBodyQuantumDot::TwoBodyQuantumDot(System* system, double alpha, double beta, double C, double omega) :
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
    m_omega = omega;
}

double TwoBodyQuantumDot::evaluate(std::vector<Particle *> particles)
{
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    for (int i=0; i<m_system->getNumberOfDimensions(); i++){
        r1 += particles[0]->getPosition()[i]*particles[0]->getPosition()[i];
        r2 += particles[1]->getPosition()[i]*particles[1]->getPosition()[i];
        r12 += (particles[0]->getPosition()[i]-particles[1]->getPosition()[i])*
                (particles[0]->getPosition()[i]-particles[1]->getPosition()[i]);
    }

    r12 = sqrt(r12);

    return m_C*exp(-m_alpha*m_omega*(r1+r2)*0.5)*exp(m_a*r12/(1+m_beta*r12));
}

double TwoBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles)
{
    double ddr = 0;
    if (m_system->getAnalyticalLaplacian()){
        //cout << "We have not yet implemented the analytically computed Laplacian." << endl;
        double r1 = 0;
        double r2 = 0;
        double r12 = 0;
        for (int i=0; i<m_system->getNumberOfDimensions(); i++){
            r1 += particles[0]->getPosition()[i]*particles[0]->getPosition()[i];
            r2 += particles[1]->getPosition()[i]*particles[1]->getPosition()[i];
            r12 += (particles[0]->getPosition()[i]-particles[1]->getPosition()[i])*
                    (particles[0]->getPosition()[i]-particles[1]->getPosition()[i]);
        }
        r12 = sqrt(r12);

        double Br12 = 1/((1+m_beta*r12)*(1+m_beta*r12));

        return m_alpha2*m_omega*m_omega*(r1+r2) - (2*m_a*Br12)*(m_alpha*m_omega*r12 - m_a*r12 - 1/r12);
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


