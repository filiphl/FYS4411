#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <iomanip>

using namespace std;

class System {
private:
    double                          m_D = 0.5;  //h2/2m. Added by us.
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 1;
    double                          m_dt = 0.005;     // Added by us.
    double                          m_derivativeStep = 1e-6;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();

    bool m_AnalyticalLaplacian = false;
    bool m_importanceSampling         = false;
    bool m_storeLocalEnergy           = false;
    bool m_storePositions             = false;
    bool m_optimizingParameters       = false;

    double prob                       = 0;
    double mynt                       = 0;
    double dx                         = 0;

public:
    bool metropolisStep             ();
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int& getNumberOfParticles()          { return m_numberOfParticles; }
    int& getNumberOfDimensions()         { return m_numberOfDimensions; }
    int& getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double& getEquilibrationFraction()   { return m_equilibrationFraction; }

    //Added by us.
    double qForce(int i, int j);
    double getStepLength() const        { return m_stepLength; }
    void openEnergyFile();
    void openPositionFile();
    void closeEnergyFile();
    void closePositionFile();
    const char* m_energyFileName   = "datafiles/localenergies.txt";
    const char* m_positionFileName = "datafiles/positions.txt";
    ofstream m_energyFile;
    ofstream m_positionFile;


    void setAnalyticalLaplacian  (bool value);
    void setImportanceSampling          (bool value);
    void setStoreLocalEnergy            (bool value);
    bool getAnalyticalLaplacian() const;
    bool getImportanceSampling() const;
    bool getStoreLocalEnergy() const;
    bool OptimizingParameters() const;
    void OptimizingParameters(bool optimizingParameters);
    bool getStorePositions() const;
    void setStorePositions(bool storePositions);
};

