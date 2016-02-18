#pragma once

class Sampler {
private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_numberOfStepsSampled    = 0;
    int     m_stepNumber              = 0;
    double  m_acceptanceRate          = 0;
    double  m_localEnergy             = 0;
    double  m_energy                  = 0;
    double  m_energySquared           = 0;
    double  m_variance                = 0;
    double  m_cumulativeEnergy        = 0;
    double  m_analyticalEnergy        = 0;
    class System* m_system = nullptr;

public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double computeAnalyticalEnergy();                     // Added by us

    // Get / Set
    double getEnergySquared() const;
    void setEnergySquared(double energySquared);

    double getVariance() const;
    void setVariance(double variance);

};
