#include <iostream>
#include <iomanip>
//#include <mpi.h>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactinsimplegaussian.h"
#include "WaveFunctions/heliumwavefunction.h"
#include "WaveFunctions/twobodyquantumdot.h"
#include "WaveFunctions/manybodyquantumdotwavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/interactingharmonicoscillator.h"
#include "Hamiltonians/heliumhamiltonian.h"
#include "Hamiltonians/twobodyquantumdothamiltonian.h"
#include "Hamiltonians/manybodyquantumdothamiltonian.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "time.h"
#include "optimizer.h"

using namespace std;

int main(int argc, char* argv[]) {



//    MPI_Init (&argc, &argv);	/* starts MPI */
//    int rank, size;
//    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
//    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
//    printf( "Hello world from process %d of %d\n", rank, size );
//    MPI_Finalize();


    int numberOfParticles   = 2;
    int numberOfDimensions  = 2;
    int numberOfSteps       = (int) 1e4;
    double omegaHO          = .5;          // Oscillator frequency.
    double omegaZ           = 1.0;
    double alpha            = 0.95455;//1.843;          // Variational parameter.
    double beta             = 0.50905;      // Variational parameter.
    double gamma            = 2.82843;
    double stepLength       = 1.3;            // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    double C                = 1.0;
    double a                = 1.0;

    System* system = new System();
    system->setInitialState                 (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setHamiltonian                  (new TwoBodyQuantumDotHamiltonian(system, omegaHO));
    //system->setWaveFunction                 (new SimpleGaussian(system, alpha));
    //system->setHamiltonian                  (new TwoBodyQuantumDotHamiltonian(system, omegaHO));
    system->setWaveFunction                 (new TwoBodyQuantumDot(system, alpha, beta, C, omegaHO, a));
    system->setEquilibrationFraction        (equilibration);
    system->setStepLength                   (stepLength);
    system->setAnalyticalLaplacian          (true);
    system->setImportanceSampling           (false);
    system->setStoreLocalEnergy             (false);
    system->setStorePositions               (false);

    //Optimizer* myOptimizer = new Optimizer(system);
    //myOptimizer->optimizeParameters();


    double t0 = clock();
    system->runMetropolisSteps              (numberOfSteps);
    double t1 = clock()-t0;


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



