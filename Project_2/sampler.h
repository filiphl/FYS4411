#pragma once

class Sampler {
private:
    class System* m_system = nullptr;

    int     m_numberOfMetropolisSteps    = 0;
    int     m_numberOfStepsSampled       = 0;
    int     m_accepted                   = 0;
    int     m_stepNumber                 = 0;
    double  m_acceptanceRate             = 0;
    double  m_localEnergy                = 0;
    double  m_energy                     = 0;
    double  m_energySquared              = 0;
    double  m_variance                   = 0;
    double  m_cumulativeEnergy           = 0;
    double  m_analyticalEnergy           = 0;

    // Optimization.
    double  m_psiDerivative              = 0;
    double  m_psiDerivativeBeta          = 0;
    double  m_cumulativePsiDeriv         = 0;
    double  m_cumulativePsiDerivBeta     = 0;
    double  m_cumulativePsiLocalProd     = 0;
    double  m_cumulativePsiLocalProdBeta = 0;
    double  m_localAlphaDeriv            = 0;
    double  m_localBetaDeriv             = 0;

public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double computeAnalyticalEnergy();                     // Added by us
    void reset();
    // Get / Set
    double getEnergySquared() const;
    void setEnergySquared(double energySquared);

    double getVariance() const;
    void setVariance(double variance);

    void setStepNumber(int stepNumber);

    double getLocalAlphaDeriv() const;
    double getLocalBetaDeriv() const;
};
