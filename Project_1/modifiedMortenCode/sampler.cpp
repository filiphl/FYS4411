#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */


    if (acceptedStep){
        m_localEnergy = m_system->getHamiltonian()->
                computeLocalEnergy(m_system->getParticles());

        m_numberOfStepsSampled++;
        m_cumulativeEnergy  += m_localEnergy;
        m_energySquared     += m_localEnergy*m_localEnergy;
    }

    if (m_system->storeLocalEnergy){
        m_system->m_outfile << setw(20) << setprecision(13) << m_localEnergy << endl;
    }

    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa[i] << endl;
    }
    cout << endl;
    cout << "  ----- Reults ----- \n" << endl;
    cout << setw(25) << left << "Numerical energy" << left << setw(25) << "Analytical energy" << endl;
    cout << setw(25) << left << m_energy           << left << setw(25) << m_analyticalEnergy  << endl<<endl;
    cout << "Variance in ergy measurements : " << m_variance << endl;
    cout << "Acceptance rate : " << setprecision(6) << m_acceptanceRate << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    m_energy         = m_cumulativeEnergy / (double)m_numberOfStepsSampled; // It is now.
    m_energySquared  = m_energySquared / (double)m_numberOfStepsSampled;
    m_variance       = m_energySquared - m_energy*m_energy;
    m_acceptanceRate = m_numberOfStepsSampled/((double)m_numberOfMetropolisSteps*(1-m_system->getEquilibrationFraction()));
}

double Sampler::computeAnalyticalEnergy()
{
    m_analyticalEnergy = m_system->getHamiltonian()->computeAnalyticalEnergy(m_system->getParticles());
}

double Sampler::getEnergySquared() const
{
    return m_energySquared;
}

void Sampler::setEnergySquared(double energySquared)
{
    m_energySquared = energySquared;
}

double Sampler::getVariance() const
{
    return m_variance;
}

void Sampler::setVariance(double variance)
{
    m_variance = variance;
}

/*
double Solver::Analytical(){
    double energy = 0;
    double h2 = 1;
    for (int p=0; p<m_nParticles; p++){
        double r_single = 0;
        for (int d=0; d<m_nDimensions; d++){
            r_single += r(p,d)*r(p,d);
        }
        energy += h2*m_alpha*(-2*m_alpha*r_single + m_nDimensions) + 0.5*m_m*m_w*m_w*r_single;
    }
    return energy;
    */

