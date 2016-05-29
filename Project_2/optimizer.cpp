#include "optimizer.h"


Optimizer::Optimizer(System *system, double alpha, double beta)
{
    m_system = system;
    m_alpha  = alpha;
    m_beta   = beta;
    system->setPrintResults(false);
}

void Optimizer::optimizeParameters()
{
    int nSteps = 1e5;
    m_system->getWaveFunction()->setAlpha(m_alpha);
    m_system->getWaveFunction()->setBeta(m_beta);

    m_system->OptimizingParameters(true);
    m_system->runMetropolisSteps(nSteps);
    m_energy = m_system->getSampler()->getEnergy();

    m_dAlpha = m_system->getSampler()->getLocalAlphaDeriv();
    m_dBeta = m_system->getSampler()->getLocalBetaDeriv();

    cout << "Alpha derivative: " << m_dAlpha  << endl;
    cout << "Beta derivative:  " << m_dBeta   << endl;
    cout << "epsilon:          " << m_epsilon << endl;

    while (true){   // fabs(m_dBeta) > m_epsilon
        m_alphaOld  = m_alpha;
        m_betaOld   = m_beta;
        m_dAlphaOld = m_dAlpha;
        m_dBetaOld  = m_dBeta;
        m_energyOld = m_energy;

        m_system->getWaveFunction()->setAlpha(m_alpha);
        m_system->getWaveFunction()->setBeta (m_beta);
        m_system->runMetropolisSteps(nSteps);
        m_energy = m_system->getSampler()->getEnergy();

        m_dAlpha = m_system->getSampler()->getLocalAlphaDeriv();
        m_dBeta  = m_system->getSampler()->getLocalBetaDeriv();
        m_system->getSampler()->setStepNumber(0);


        m_alpha -= m_steplength*m_dAlpha;
        m_beta  -= m_steplengthBeta*m_dBeta;
        if (m_dBeta*m_dBeta + m_dAlpha*m_dAlpha > m_dBetaOld*m_dBetaOld + m_dAlphaOld*m_dAlphaOld){   //

/*
            if (m_dAlpha*m_dAlpha > m_dAlphaOld*m_dAlphaOld){
                m_steplength *= 0.9;
            }
            if (m_dBeta*m_dBeta > m_dBetaOld*m_dBetaOld){
                m_steplengthBeta *= 0.9;
            }
*/
            m_steplength *= 0.7;
            m_steplengthBeta *= 0.7;

            m_dAlpha = m_dAlphaOld;
            m_dBeta  = m_dBetaOld;
            m_alpha  = m_alphaOld;
            m_beta   = m_betaOld;
            nSteps = (int)nSteps*1.0;
            if (m_steplength < 1e-6 || m_steplengthBeta < 1e-6) {
                m_system->runMetropolisSteps(nSteps);
                m_energy = m_system->getSampler()->getEnergy();
                cout << "alpha                        : " << m_alphaOld      << endl;
                cout << "beta                         : " << m_betaOld       << endl;
                cout << "Derivative step length alpha : " << m_steplength << endl;
                cout << "Derivative step length beta  : " << m_steplengthBeta << endl;
                cout << "Beta derivative              : " << m_dBetaOld      << endl;
                cout << "Energy                       : " << m_energyOld<<endl<<endl;
                return;
            }
        }
        cout << "alpha                        : " << setprecision(10) << m_alpha << endl;
        cout << "beta                         : " << setprecision(10) << m_beta  << endl;
        cout << "Derivative step length       : " << m_steplength << endl;
        cout << "Derivative step length beta  : " << m_steplengthBeta << endl;
        cout << "Alpha derivative             : " << m_dAlpha << endl;
        cout << "Beta derivative              : " << m_dBeta  << endl;
        cout << "Energy                       : " << m_system->getSampler()->getEnergy()<<endl<<endl;
        m_system->getSampler()->reset();
    }
}


double Optimizer::absoluteValue(double value){ return sqrt(value*value); }


void Optimizer::setAlpha     (double alpha)     { m_alpha           = alpha;     }
void Optimizer::setBeta      (double beta)      { m_beta            = beta;      }
void Optimizer::setEpsilon   (double epsilon)   { m_epsilon         = epsilon;   }
//void Optimizer::setSteplength(double steplength){ m_steplengthAlpha = steplength;}

double Optimizer::getAlpha()   const { return m_alpha;           }
double Optimizer::getBeta()    const { return m_beta;            }
double Optimizer::epsilon()    const { return m_epsilon;         }
//double Optimizer::steplength() const { return m_steplengthAlpha; }





