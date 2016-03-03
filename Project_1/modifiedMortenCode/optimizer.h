#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "system.h"

class Optimizer
{
private:
    double m_steplength     = 0.1;
    double m_epsilon        = 1e-4;
    double m_alpha          = 0.4;
    double m_alphaOld       = 0.4;
    double m_derivative     = 0;
    double m_derivativeOld  = 0;
    int    m_numberOfSteps  = 0;
    System* m_system        = nullptr;
public:
    Optimizer(System* system);
    void optimizeParameters();
    double absoluteValue(double value);
    double getAlpha() const;
    void setAlpha(double alpha);
    double epsilon() const;
    void setEpsilon(double epsilon);
    double steplength() const;
    void setSteplength(double steplength);
};

#endif // OPTIMIZER_H
