#pragma once
#include <vector>
#include "system.h"
#include <armadillo>

class WaveFunction {


protected:
    int m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    double m_alpha  = 0;
    double m_alpha2 = 0;
    double m_beta   = 0;
    double m_beta2  = 0;
    int m_newlyMoved = 0;

public:
    arma::mat    m_distances;
    WaveFunction(){}
    WaveFunction(class System* system);
    virtual double evaluate(std::vector<class Particle*> particles)         = 0;
    virtual double computeLaplacian(std::vector<class Particle*> particles) = 0;
    virtual double computeGradient(std::vector<class Particle*> particles, int particle, int dimension) = 0;
    virtual double computeRatio(std::vector<class Particle*> particles, int i, int j, double change)    = 0;
    virtual void   updateSlater(int i) = 0;

    double sumOfArguments = 0;

    double getAlpha() const;
    void setAlpha(double alpha);
    double getAlpha2() const;
    void setAlpha2(double alpha2);
    int getNumberOfParameters()         { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    double getBeta() const;
    void setBeta(double beta);
    double getBeta2() const;
    void setBeta2(double beta2);
    void setNewlyMoved(int i) { m_newlyMoved = i; }
};

