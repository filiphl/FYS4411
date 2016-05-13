#include "twobodyquantumdot.h"

TwoBodyQuantumDot::TwoBodyQuantumDot(System* system, double alpha, double beta, double C, double omega, double a) :
    WaveFunction(system)
{
    if (system->getNumberOfParticles() != 2){
        cout << "The two body quantum dot contains 2 electrons (particles)." << endl;
        cout << "Number of particles given: "<< system->getNumberOfParticles() << endl;
        exit(0);
    }
    assert(alpha >= 0);
    m_numberOfParameters = 5;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(C);
    m_parameters.push_back(omega);
    m_parameters.push_back(a);
    setAlpha(alpha);
    setBeta(beta);
    setAlpha(alpha);
    m_a = a;
    m_C = C;
    m_omega = omega;
    name = "Two body quantum dot";
}

double TwoBodyQuantumDot::evaluate(std::vector<Particle *> particles)
{
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    for (int i=0; i<m_system->getNumberOfDimensions(); i++){
        r1  +=  particles[0]->getNewPosition()[i]*particles[0]->getNewPosition()[i];
        r2  +=  particles[1]->getNewPosition()[i]*particles[1]->getNewPosition()[i];
        r12 += (particles[0]->getNewPosition()[i]-particles[1]->getNewPosition()[i])*
               (particles[0]->getNewPosition()[i]-particles[1]->getNewPosition()[i]);
    }


    r12 = sqrt(r12);

    psiAlpha = -m_alpha*m_omega*(r1+r2)*0.5;
    psiBeta  = -m_a*r12*r12/((1+m_beta*r12)*(1+m_beta*r12));
    cout << exp(m_a*r12/(1+m_beta*r12))<<endl;
    return m_C*exp(-m_alpha*m_omega*(r1+r2)*0.5)*exp(m_a*r12/(1+m_beta*r12));
}

double TwoBodyQuantumDot::computeLaplacian(std::vector<Particle *> particles)
{
    double ddr = 0;
    if (m_system->getAnalyticalLaplacian()){
//        cout << "Analytical"<<endl;
        // return the laplacian/psi
        double r1 = 0;
        double r2 = 0;
        double r12 = 0;
        for (int i=0; i<m_system->getNumberOfDimensions(); i++){
            r1  +=  particles[0]->getNewPosition()[i]*particles[0]->getNewPosition()[i];
            r2  +=  particles[1]->getNewPosition()[i]*particles[1]->getNewPosition()[i];
            r12 += (particles[0]->getNewPosition()[i]-particles[1]->getNewPosition()[i])*
                   (particles[0]->getNewPosition()[i]-particles[1]->getNewPosition()[i]);
        }
        r12 = sqrt(r12);

        double Br12 = (1+m_beta*r12);

        return - 2*m_system->getNumberOfDimensions()*m_alpha*m_omega
               - (4*m_a*m_beta)/(Br12*Br12*Br12)
               + 2*(m_system->getNumberOfDimensions()-1)*m_a/(Br12*Br12*r12)
               + m_alpha2*m_omega*m_omega*(r1+r2)
               + (2*m_a*m_a)/(Br12*Br12*Br12*Br12)
               - (2*m_alpha*m_omega*m_a*r12)/(Br12*Br12);
    }

    else {
//        cout << "Numerical"<<endl;
        m_derivativeStepLength = 1e-5;
        const double psi = evaluate( particles );
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                particles[i]->adjustNewPosition( m_derivativeStepLength, j );      // +
                const double psiPlus  =   evaluate( particles );
                particles[i]->adjustNewPosition( -2 * m_derivativeStepLength, j ); // -
                const double psiMinus =   evaluate( particles );
                particles[i]->adjustNewPosition( m_derivativeStepLength, j );      // reset
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr / (m_derivativeStepLength * m_derivativeStepLength);
    }
    return ddr;
}


double TwoBodyQuantumDot::computeGradient(std::vector<Particle *> particles, int particle , int dimension)
{
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    for (int i=0; i<m_system->getNumberOfDimensions(); i++){
        r1  +=  particles[0]->getNewPosition()[i]*particles[0]->getNewPosition()[i];
        r2  +=  particles[1]->getNewPosition()[i]*particles[1]->getNewPosition()[i];
        r12 += (particles[0]->getNewPosition()[i]-particles[1]->getNewPosition()[i])*
               (particles[0]->getNewPosition()[i]-particles[1]->getNewPosition()[i]);
    }
    r12 = sqrt(r12);
    double Br12 = (1+m_beta*r12)*(1+m_beta*r12);

    return -m_alpha*m_omega*particles[particle]->getNewPosition()[dimension]
            + (m_a*(1-2*(particle%2==0)*(particles[0]->getNewPosition()[dimension]-particles[1]->getNewPosition()[dimension])))
            /(r12*Br12);
}

double TwoBodyQuantumDot::computeRatio(std::vector<class Particle*> particles, int i, int j, double change)
{
    double oldPsi = evaluate(particles);
    particles[i]->adjustNewPosition(change, j);
    double newPsi = evaluate(particles);
    //particles[i]->adjustOldPosition(-change, j);    // This was the problem Sean. The twoBodyQuantumDot got moved twice, if accepted. See mynt<prob in system.cpp
    return newPsi/oldPsi;
}


double TwoBodyQuantumDot::computeDerivativeOfAlpha(){
    double sum = 0;
    cout << "YES"<<endl;
    for (Particle* particle : m_system->getParticles()){
        for (int i=0; i< m_system->getNumberOfDimensions(); i++){
            sum -= particle->getNewPosition()[i]*particle->getNewPosition()[i];

        }
    }
    return sum;
}

std::vector<double> TwoBodyQuantumDot::getParameters()
{

}

void TwoBodyQuantumDot::updateSlater(int i)
{

}

void TwoBodyQuantumDot::printParameters()
{

}


