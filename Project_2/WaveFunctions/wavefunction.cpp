#include "wavefunction.h"
#include <iostream>     // Added by us

double WaveFunction::getAlpha() const
{
    return m_alpha;

}

void WaveFunction::setAlpha(double alpha)
{
    m_alpha  = alpha;
    m_alpha2 = alpha*alpha;
}

double WaveFunction::getAlpha2() const
{
    return m_alpha2;
}

void WaveFunction::setAlpha2(double alpha2)
{
    m_alpha2 = alpha2;
}

double WaveFunction::getBeta2() const
{
    return m_beta2;
}

void WaveFunction::setBeta2(double beta2)
{
    m_beta2 = beta2;
}

double WaveFunction::getBeta() const
{
    return m_beta;
}

void WaveFunction::setBeta(double beta)
{
    m_beta  = beta;
    m_beta2 = beta*beta;
}

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}




