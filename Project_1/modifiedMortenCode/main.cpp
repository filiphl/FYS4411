#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "time.h"

using namespace std;

int main() {
    int numberOfParticles   = 10;           // This is the number of particles. P.A.R.T.I.C.L.E.S.
    int numberOfDimensions  = 3;
    int numberOfSteps       = (int) 1e4;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used
                                            // for equilibration.

    System* system = new System();
    system->setHamiltonian(new HarmonicOscillator(system, omega));
    system->setWaveFunction                 (new SimpleGaussian(system, alpha));
    system->setInitialState                 (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction        (equilibration);
    system->setStepLength                   (stepLength);
    system->setAnalyticalDoubleDerivative   (true);
    system->setStoreLocalEnergy             (false);
    system->setImportanceSampling           (true);
    double t0 = clock();
    system->runMetropolisSteps              (numberOfSteps);
    double t1 = clock()-t0;
    cout << t1*1e-6 << endl;
    return 0;
}
