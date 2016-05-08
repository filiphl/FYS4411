#include "optimizer.h"


Optimizer::Optimizer(System *system)
{
    m_system = system;
}

void Optimizer::optimizeParameters()
{
    int nSteps = 1e3;
    m_system->getWaveFunction()->setAlpha(m_alpha);
    m_system->OptimizingParameters(true);
    m_system->runMetropolisSteps(nSteps);

    m_dAlpha     = m_system->getSampler()->getLocalAlphaDeriv();
    m_dBeta = m_system->getSampler()->getLocalBetaDeriv();
    cout << "derivative alpha: "<< m_dAlpha    << endl;
    cout << "derivative beta: " << m_dBeta << endl;
    cout << "epsilon: "<< m_epsilon<< endl;

    while (true){ //(fabs(m_dAlpha) > m_epsilon){
        m_alphaOld  = m_alpha;
        m_dAlphaOld = m_dAlpha;
        m_betaOld   = m_beta;
        m_dBetaOld  = m_dBeta;

        m_alpha -= m_steplengthAlpha*m_dAlpha;
        m_beta  -= m_steplengthBeta*m_dBeta;
        m_system->getWaveFunction()->setAlpha(m_alpha);
        m_system->getWaveFunction()->setBeta(m_beta);
        m_system->runMetropolisSteps(nSteps);
        m_dAlpha = m_system->getSampler()->getLocalAlphaDeriv();
        m_dBeta  = m_system->getSampler()->getLocalBetaDeriv();
        m_system->getSampler()->setStepNumber(0);
//        cout << m_dAlpha*m_dAlpha + m_dBeta*m_dBeta<<"    "<<m_dAlphaOld*m_dAlphaOld + m_dBetaOld*m_dBetaOld<<endl;
        if (m_dAlpha*m_dAlpha + m_dBeta*m_dBeta > m_dAlphaOld*m_dAlphaOld + m_dBetaOld*m_dBetaOld){
            if (fabs(m_dAlpha) > fabs(m_dAlphaOld)){ m_steplengthAlpha /= 1.2; }
            if (fabs(m_dBeta)  > fabs(m_dBetaOld)) { m_steplengthBeta  /= 1.2; }
            //cout << "m_dBeta = "<< fabs(m_dBeta) << "   m_dBetaOld = "<< m_dBetaOld<<endl;
            //cout << "m_dAlpha = "<< fabs(m_dAlpha) << "   m_dAlphaOld = "<< m_dAlphaOld<<endl;
            nSteps        = (int)nSteps*1;
            if (m_steplengthAlpha < 1e-6) {
                cout << "Energy                        " << m_system->getSampler()->getEnergy()<<endl;
                cout << "alpha                         " << m_alpha           << endl;
                cout << "beta                          " << m_beta            << endl;
                cout << "Derivative step length alpha: " << m_steplengthAlpha << endl;
                cout << "Derivative step length beta:  " << m_steplengthAlpha << endl;
                cout << "Alpha derivative:             " << m_dAlpha          << endl;
                cout << "Beta  derivative:             " << m_dBeta           << endl;
                cout << "\n-------------------\n";
                return;
            }
            m_dAlpha      = m_dAlphaOld;
            m_dBeta       = m_dBetaOld;
            m_alpha       = m_alphaOld;
            m_beta        = m_betaOld;
        }
        else{
            cout << "Energy                        " << m_system->getSampler()->getEnergy()<<endl;
            cout << "alpha                         " << m_alpha           << endl;
            cout << "beta                          " << m_beta            << endl;
            cout << "Derivative step length alpha: " << m_steplengthAlpha << endl;
            cout << "Derivative step length beta:  " << m_steplengthAlpha << endl;
            cout << "Alpha derivative:             " << m_dAlpha          << endl;
            cout << "Beta  derivative:             " << m_dBeta           << endl;
            cout << "\n-------------------\n";

        }
        m_system->getSampler()->reset();
    }
}

double Optimizer::absoluteValue(double value)
{
    return sqrt(value*value);
}

double Optimizer::getAlpha() const
{
    return m_alpha;
}

void Optimizer::setAlpha(double alpha)
{
    m_alpha = alpha;
}

double Optimizer::epsilon() const
{
    return m_epsilon;
}

void Optimizer::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

double Optimizer::steplength() const
{
    return m_steplengthAlpha;
}

void Optimizer::setSteplength(double steplength)
{
    m_steplengthAlpha = steplength;
}
