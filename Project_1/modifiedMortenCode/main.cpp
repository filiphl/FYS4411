#include <iostream>
#include <iomanip>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactinsimplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/interactingharmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "time.h"
#include "optimizer.h"

using namespace std;

int main() {



    int numberOfParticles   = 10;            // This is the number of particles. P.A.R.T.I.C.L.E.S.
    int numberOfDimensions  = 3;
    int numberOfSteps       = (int) 1e5;
    double omegaHO          = 1.0;          // Oscillator frequency.
    double omegaZ           = 1.0;
    double alpha            = 0.49746;          // Variational parameter.
    double beta             = 2.82843;      // Variational parameter.
    double gamma            = 2.82843;
    double stepLength       = 1;            // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.

    System* system = new System();
    system->setHamiltonian                  (new InteractingHarmonicOscillator(system, omegaHO, omegaZ, gamma));//(new HarmonicOscillator(system, omegaHO));
    system->setWaveFunction                 (new InteractinSimpleGaussian(system, alpha, beta));//(new SimpleGaussian(system, alpha));
    //system->setHamiltonian                  (new HarmonicOscillator(system, omegaHO));
    //system->setWaveFunction                 (new SimpleGaussian(system, alpha));
    system->setInitialState                 (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction        (equilibration);
    system->setStepLength                   (stepLength);
    system->setAnalyticalDoubleDerivative   (true);
    system->setImportanceSampling           (false);
    system->setStoreLocalEnergy             (false);
    system->setStorePositions               (true);

    Optimizer* myOptimizer = new Optimizer(system);
    //myOptimizer->optimizeParameters();


    double t0 = clock();
    system->runMetropolisSteps              (numberOfSteps);
    double t1 = clock()-t0;
    cout << setprecision(8) << t1*1e-6 << endl;
    return 0;
}


/*

class Optimizer {
public:
    Optimizer(System* system);
    double runOptimize();


private:
    System* m_system = nullptr;

};



int main() {
    System mySystem = new System();
    mySystem.setWaveFunction(new InteractinSimpleGaussian());
    ...

    Optimizer myOptimizer = new Optimizer(mySystem);
    myOptimizer.runOptimize();


}


double Optimizer::runOptimize() {
    while derivative > eps
            m_system->runMetropolisSteps(10^6);
            m_system->getSampler()->hva var alpha derivative?
            m_system->adjustAlpha(.,.)
    return alpha;
}
  */



