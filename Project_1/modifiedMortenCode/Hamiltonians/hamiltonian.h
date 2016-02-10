#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    virtual double computeAnalyticalEnergy(std::vector<class Particle*> particles) = 0;
    double computeKineticEnergy(std::vector<Particle *> particles);
protected:
    class System* m_system = nullptr;
    class WaveFunction* m_waveFunction = nullptr;
};

