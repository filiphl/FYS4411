#include <iostream>
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

using namespace std;

int main() {
    int numberOfParticles   = 5;            // This is the number of particles. P.A.R.T.I.C.L.E.S.
    int numberOfDimensions  = 3;
    int numberOfSteps       = (int) 1e3;
    double omegaHO          = 1.0;          // Oscillator frequency.
    double omegaZ           = 1.0;
    double alpha            = 0.5;          // Variational parameter.
    double beta             = 1;            // Variational parameter.
    double gamma            = 2.82843;
    double stepLength       = 1;            // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.



    System* system = new System();
    system->setHamiltonian                  (new InteractingHarmonicOscillator(system, omegaHO, omegaZ, gamma));//(new HarmonicOscillator(system, omegaHO));
    system->setWaveFunction                 (new InteractinSimpleGaussian(system, alpha, beta));//(new SimpleGaussian(system, alpha));
    system->setInitialState                 (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction        (equilibration);
    system->setStepLength                   (stepLength);
    system->setAnalyticalDoubleDerivative   (false);
    system->setStoreLocalEnergy             (true);
    system->setImportanceSampling           (false);
    double t0 = clock();
    system->runMetropolisSteps              (numberOfSteps);
    double t1 = clock()-t0;
    cout << t1*1e-6 << endl;
    return 0;
}
