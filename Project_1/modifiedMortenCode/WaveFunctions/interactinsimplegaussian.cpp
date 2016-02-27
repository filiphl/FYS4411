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

        for (int k=0; k< m_system->getNumberOfParticles(); k++){
            double term1 = 0;
            double term2 = 0;
            double term3 = 0;
            double term4 = 0;
            double term5 = 0;
            Particle* particleK = m_system->getParticles()[k];

            // Term1
            for (int d=0; d<m_system->getNumberOfDimensions(); d++){
                if (d<2){
                    term1 += particleK->getPosition()[d]*particleK->getPosition()[d];
                }
                else{
                    term1 += particleK->getPosition()[d]*particleK->getPosition()[d]*m_beta2;
                }
            }

            // Term 2
            for (int i=0; i<m_system->getNumberOfParticles(); i++){
                if (i!=k){
                    double temp2 = 0;
                    Particle* particleI = m_system->getParticles()[i];

                    for (int d=0; d<m_system->getNumberOfDimensions(); d++){
                        if (d<2){
                            temp2 += particleK->getPosition()[d] * ( particleI->getPosition()[d] - particleK->getPosition()[d] );
                        }
                        else{
                            temp2 += particleK->getPosition()[d] * ( particleI->getPosition()[d] - particleK->getPosition()[d] ) * m_beta;
                        }
                    }
                    temp2 *= uOverR(particleI, particleK);
                    term2 += temp2;

                    // Remember to devide by r_ki and multiply by u'

                    //Term3
                    for (int j=0; j<m_system->getNumberOfParticles(); j++){
                        if (j!=k){
                            Particle* particleJ = m_system->getParticles()[j];
                            double temp3 = 0;
                            for (int d=0; d<m_system->getNumberOfDimensions(); d++){
                                temp3 += (particleK->getPosition()[d]-particleI->getPosition()[d]) *
                                         (particleK->getPosition()[d]-particleJ->getPosition()[d]);
                            }
                            temp3 *= uOverR(particleK, particleJ);
                            temp3 *= uOverR(particleI, particleK);
                            term3 += temp3;
                        }
                    }



                    // Term4

                    term4 += 2*uOverR(particleI, particleK);         // (d-1)


                    // Term5

                    term5 += u2OverR2(particleI, particleK);
                }
            }



            term1 *= 4*m_alpha2;
            term1 -= 2*m_alpha*(2+m_beta);
            term2 *= -4*m_alpha;
            ddr += term1 + term2 + term3 + term4 + term5;

            cout << "  term1 " << setw(10) << setprecision(4) << left <<term1;
            cout << "  term2 " << setw(10) << setprecision(4) << left <<term2;
            cout << "  term3 " << setw(10) << setprecision(4) << left <<term3;
            cout << "  term4 " << setw(10) << setprecision(4) << left <<term4;
            cout << "  term5 " << setw(10) << setprecision(4) << left <<term5 << endl;

        }
        return ddr;
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




double InteractinSimpleGaussian::uOverR(Particle* particle1, Particle* particle2){
    double r = 0;
    double u = 0;
    for (int i=0; i<m_system->getNumberOfDimensions(); i++){
        r += (particle1->getPosition()[i] - particle2->getPosition()[i]) *
                (particle1->getPosition()[i] - particle2->getPosition()[i]);
    }
    r = sqrt(r);
    u = m_a / ( r*r - m_a*r );

    return u/r;
}



double InteractinSimpleGaussian::u2OverR2(Particle *particle1, Particle *particle2)
{
    double r = 0;
    double u = 0;
    for (int i=0; i<m_system->getNumberOfDimensions(); i++){
        r += (particle1->getPosition()[i] - particle2->getPosition()[i]) *
                (particle1->getPosition()[i] - particle2->getPosition()[i]);
    }
    r = sqrt(r);
    return -m_a*(2*r-m_a)/(r*r-m_a*r);
}







//  Get / Set


double InteractinSimpleGaussian::alpha() const
{
    return m_alpha;
}

void InteractinSimpleGaussian::setalpha(double alpha)
{
    m_alpha = alpha;
    m_alpha2 = alpha*alpha;
}

double InteractinSimpleGaussian::beta() const
{
    return m_beta;
}

void InteractinSimpleGaussian::setBeta(double beta)
{
    m_beta = beta;
    m_beta2 = beta*beta;
}
