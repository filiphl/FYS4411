#pragma once
#include <vector>


class WaveFunction {


protected:
    int m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    double m_alpha  = 0;
    double m_alpha2 = 0;
    double m_beta   = 0;
    double m_beta2  = 0;

public:
    WaveFunction(){};
    WaveFunction(class System* system);
    virtual double evaluate(std::vector<class Particle*> particles)         = 0;
    virtual double computeLaplacian(std::vector<class Particle*> particles) = 0;
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
};

