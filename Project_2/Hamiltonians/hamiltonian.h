#pragma once
#include "../system.h"

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    virtual double computeAnalyticalEnergy(std::vector<class Particle*> particles) = 0;
    double computeKineticEnergy(std::vector<Particle *> particles);

    double getKineticEnergy() const;

    double getPotentialEnergy() const;

protected:
    class System*       m_system       = nullptr;
    class WaveFunction* m_waveFunction = nullptr;
    double m_kineticEnergy   = 0;
    double m_potentialEnergy = 0;
};

