#include "interactinsimplegaussian.h"


InteractinSimpleGaussian::InteractinSimpleGaussian(System *system, double alpha, double beta, double gamma) :
    WaveFunction(system)
{
    assert(alpha >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(gamma);
    setalpha(alpha);
    setBeta(beta);
    setGamma(gamma);
}

double InteractinSimpleGaussian::evaluate(std::vector<Particle *> particles)
{
    double argument = 0;
    for (int i=0; i< m_system->getNumberOfParticles(); i++){
        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
            if (j<2){
                argument -= particles[i]->getPosition()[j] * particles[i]->getPosition()[j];
            }
            else{
                argument -= m_beta*particles[i]->getPosition()[j] * particles[i]->getPosition()[j];
            }
        }
    }
    g = exp(m_alpha * argument);

    double f = 1;
    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        for (int j=i+1; j<m_system->getNumberOfParticles(); j++){
            for (int k=0; k<m_system->getNumberOfDimensions(); k++){
                dr2 = (particles[j]-particles[i]) * (particles[j]-particles[i]);
            }
            absdr = sqrt(dr2);
            if (absdr <= m_a){ return 0; }
            else{ f *= 1-a/absdr; }
        }
    }

    return g*f;
}

void InteractinSimpleGaussian::computeDoubleDerivative(std::vector<Particle *> particles)
{

}











//  Get / Set


double InteractinSimpleGaussian::alpha() const
{
    return m_alpha;
}

void InteractinSimpleGaussian::setalpha(double alpha)
{
    m_alpha = alpha;
}

double InteractinSimpleGaussian::beta() const
{
    return m_beta;
}

void InteractinSimpleGaussian::setBeta(double beta)
{
    m_beta = beta;
}

double InteractinSimpleGaussian::gamma() const
{
    return m_gamma;
}

void InteractinSimpleGaussian::setGamma(double gamma)
{
    m_gamma = gamma;
}
