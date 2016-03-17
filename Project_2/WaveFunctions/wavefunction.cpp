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

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}




