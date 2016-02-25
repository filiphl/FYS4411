#include "interactinsimplegaussian.h"


InteractinSimpleGaussian::InteractinSimpleGaussian(System *system, double alpha, double beta) :
    WaveFunction(system)
{
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    cout << m_numberOfParameters<<endl;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    setalpha(alpha);
    setBeta(beta);
}

double InteractinSimpleGaussian::evaluate(std::vector<Particle *> particles)
{
    double argument = 0;
    for (int i=0; i< m_system->getNumberOfParticles(); i++){
        for (int k=0; k<m_system->getNumberOfDimensions(); k++){
            if (k<2){
                argument -= particles[i]->getPosition()[k] * particles[i]->getPosition()[k];
            }
            else{
                argument -= m_beta*particles[i]->getPosition()[k] * particles[i]->getPosition()[k];
            }
        }
    }
    double g = exp(m_alpha * argument);

    double   f = 1;
    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        for (int j=i+1; j<m_system->getNumberOfParticles(); j++){
            double dr2 = 0;
            for (int k=0; k<m_system->getNumberOfDimensions(); k++){
                dr2 += (particles[i]->getPosition()[k] - particles[j]->getPosition()[k]) *
                       (particles[i]->getPosition()[k] - particles[j]->getPosition()[k]) ;
            }
            absdr = sqrt(dr2);
            if (absdr <= m_a){ return 0; }
            else{ f *= 1-m_a/absdr; }
        }
    }
    return g*f;
}

double InteractinSimpleGaussian::computeDoubleDerivative(std::vector<Particle *> particles)
{
    double ddr = 0;

    if (m_system->getAnalyticalDoublederivative()){
        cout << "We have not yet implemented the analytical solution to the double derivative, though we should."<< endl;
        exit(1);
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
    return ddr;
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
