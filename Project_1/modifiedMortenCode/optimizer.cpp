#include "optimizer.h"


Optimizer::Optimizer(System *system)
{
    m_system = system;
}

void Optimizer::optimizeParameters()
{
    int nSteps = 1e4;
    m_system->getWaveFunction()->setAlpha(m_alpha);
    m_system->OptimizingParameters(true);
    m_system->runMetropolisSteps(nSteps);
    m_derivative = m_system->getSampler()->getLocalAlphaDeriv();
    cout << "derivative: "<<m_derivative << endl;
    cout << "epsilon: "<< m_epsilon<< endl;

    while (fabs(m_derivative) > m_epsilon){
        m_derivativeOld = m_derivative;
        m_alphaOld = m_alpha;
        m_alpha -= m_steplength*m_derivative;
        m_system->getWaveFunction()->setAlpha(m_alpha);
        m_system->runMetropolisSteps(nSteps);
        m_derivative = m_system->getSampler()->getLocalAlphaDeriv();
        m_system->getSampler()->setStepNumber(0);

        if (fabs(m_derivative) > fabs(m_derivativeOld)){
            m_derivative = m_derivativeOld;
            m_steplength /= 2.;
            m_alpha = m_alphaOld;
            nSteps = (int)nSteps*1.1;
            if (m_steplength < 1e-7) {
                cout << "alpha " << m_alpha << endl;
                cout << "Derivative step length: " << m_steplength << endl;
                cout << "Alpha derivative: " << m_derivative << endl;
                return;
            }
        }
        cout << "alpha " << setprecision(10)<< m_alpha << endl;
        cout << "Derivative step length: " << m_steplength << endl;
        cout << "Alpha derivative: " << m_derivative << endl;
    }
    m_system->getSampler()->reset();
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
    return m_steplength;
}

void Optimizer::setSteplength(double steplength)
{
    m_steplength = steplength;
}
