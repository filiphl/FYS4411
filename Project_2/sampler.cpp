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

struct Pos {
    double x,y;
} pos;

double Sampler::getLocalAlphaDeriv() const
{
    return m_localAlphaDeriv;
}

double Sampler::getLocalBetaDeriv() const
{
    return m_localBetaDeriv;
}

void Sampler::setStepNumber(int stepNumber)
{
    m_stepNumber = stepNumber;
}


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
        m_cumulativeEnergy          = 0;
        m_cumulativePsiDeriv        = 0;
        m_cumulativePsiLocalProd    = 0;
        m_accepted                  = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */



    if (acceptedStep) {
        m_accepted++;
    }

    m_localEnergy = m_system->getHamiltonian()->
            computeLocalEnergy(m_system->getParticles());

    if (m_system->OptimizingParameters()){
        m_psiDerivative     = m_system->getWaveFunction()->psiAlpha;
        m_psiDerivativeBeta = m_system->getWaveFunction()->psiBeta;
        m_cumulativePsiDeriv         += m_psiDerivative;
        m_cumulativePsiLocalProd     += m_localEnergy*m_psiDerivative;
        m_cumulativePsiDerivBeta     += m_psiDerivativeBeta;
        m_cumulativePsiLocalProdBeta += m_localEnergy*m_psiDerivativeBeta;
    }

    m_numberOfStepsSampled++;
    m_cumulativeEnergy          += m_localEnergy;
    m_cumulativeKineticEnergy   += m_system->getHamiltonian()->getKineticEnergy();
    m_cumulativePotentialEnergy += m_system->getHamiltonian()->getPotentialEnergy();
    m_energySquared             += m_localEnergy*m_localEnergy;


    if (m_system->m_energyFile.is_open()){
        m_system->m_energyFile.write(reinterpret_cast<const char*>(&m_localEnergy), sizeof(double));
    }

    if (m_system->m_oldPositionFile.is_open()){


        for (int i=0; i<m_system->getNumberOfParticles(); i++){

                pos.x = m_system->getParticles()[i]->getOldPosition()[0];
                pos.y = m_system->getParticles()[i]->getOldPosition()[1];

                m_system->m_oldPositionFile.write(reinterpret_cast<char*>(&pos), sizeof(Pos));  // Pos is a struct.
        }
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
    cout << " Name : "<<m_system->getWaveFunction()->getName()<<endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    m_system->getWaveFunction()->printParameters();
    cout << endl;
    cout << "  ----- Reults -----" << endl;
    cout << " Energy          : "  << setw(25) << setprecision(5) << left << m_energy <<endl;
    cout << " Variance        : "  << m_variance << endl;
    cout << " Acceptance rate : "  << setprecision(6) << m_acceptanceRate << endl;
    cout << " Kinetic         : "  << setw(25) << setprecision(5) << left << m_kineticEnergy   << endl;
    cout << " Potential       : "  << setw(25) << setprecision(5) << left << m_potentialEnergy << endl;
    cout << endl;

    //cout << np<<"& " << nd <<"& "<< m_analyticalEnergy<<"& " << m_energy<<"& " << m_variance<<"& ";
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    m_energy          = m_cumulativeEnergy                    / (double)m_numberOfStepsSampled;
    m_kineticEnergy   = m_cumulativeKineticEnergy             / (double)m_numberOfStepsSampled;
    m_potentialEnergy = m_cumulativePotentialEnergy           / (double)m_numberOfStepsSampled;
    m_energySquared   = m_energySquared                       / (double)m_numberOfStepsSampled;
    m_acceptanceRate  = m_accepted                            / (double)m_numberOfStepsSampled;
    m_variance        = (m_energySquared - m_energy*m_energy) / (double)m_numberOfStepsSampled;

    if (m_system->OptimizingParameters()){
        m_localAlphaDeriv = 2*(m_cumulativePsiLocalProd - m_cumulativePsiDeriv*m_energy)/(double)m_numberOfStepsSampled;
        m_localBetaDeriv  = 2*(m_cumulativePsiLocalProdBeta - m_cumulativePsiDerivBeta*m_energy)/(double)m_numberOfStepsSampled;
    }
}

double Sampler::computeAnalyticalEnergy()
{
    m_analyticalEnergy = m_system->getHamiltonian()->computeAnalyticalEnergy(m_system->getParticles());
}

void Sampler::reset()
{
    m_numberOfStepsSampled       = 0;
    m_accepted                   = 0;
    m_stepNumber                 = 0;
    m_acceptanceRate             = 0;
    m_localEnergy                = 0;
    m_psiDerivative              = 0;
    m_psiDerivativeBeta          = 0;
    m_energy                     = 0;
    m_energySquared              = 0;
    m_variance                   = 0;
    m_cumulativeEnergy           = 0;
    m_cumulativePsiDeriv         = 0;
    m_cumulativePsiDerivBeta     = 0;
    m_cumulativePsiLocalProd     = 0;
    m_cumulativePsiLocalProdBeta = 0;
    m_analyticalEnergy           = 0;
    m_localAlphaDeriv            = 0;
    m_localBetaDeriv             = 0;
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

