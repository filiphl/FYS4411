#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "system.h"

class Optimizer
{
private:
    System* m_system           = nullptr;

    double m_energy            = 0;
    double m_energyOld         = 0;
    double m_steplength        = 0.5;
    double m_steplengthBeta    = 0.5;
    double m_epsilon           = 0;
    double m_alpha             = 1;
    double m_alphaOld          = 0.0;
    double m_beta              = 0;
    double m_betaOld           = 0.0;
    double m_dAlpha            = 0;
    double m_dAlphaOld         = 0;
    double m_dBeta             = 0;
    double m_dBetaOld          = 0;
    int    m_numberOfSteps     = 0;

public:
    Optimizer(System* system, double alpha, double beta);
    void optimizeParameters();
    double absoluteValue(double value);
    double getAlpha() const;
    void setAlpha(double alpha);
    double epsilon() const;
    void setEpsilon(double epsilon);
    double steplength() const;
    void setSteplength(double steplength);
    double getBeta() const;
    void setBeta(double beta);
};

#endif // OPTIMIZER_H
