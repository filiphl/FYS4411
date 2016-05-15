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
    int numberOfSteps       = (int) 1e6;
    double omegaHO          = 1.;           // Oscillator frequency.
    double omegaZ           = 1.0;
    double alpha            = 0.94295;      // Variational parameter.
    double beta             = 0.46034;      // Variational parameter.
    double gamma            = 2.82843;
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.1;          // Fraction steps used for equilibration.
    double C                = 1.0;
    double a                = 1.0;

    System* system = new System();
    system->setInitialState                 (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    if (0){
        system->setWaveFunction                 (new ManyBodyQuantumDotWaveFunction(system, alpha, omegaHO, a, beta));
        system->setHamiltonian                  (new ManyBodyQuantumDotHamiltonian (system, omegaHO));
    }
    else{
        system->setWaveFunction                 (new TwoBodyQuantumDot(system, alpha, beta,C, omegaHO, a));
        system->setHamiltonian                  (new TwoBodyQuantumDotHamiltonian(system, omegaHO));
    }

    system->setEquilibrationFraction        (equilibration);
    system->setStepLength                   (stepLength);
    system->setAnalyticalLaplacian          (true);
    system->setImportanceSampling           (true);
    system->setStoreLocalEnergy             (false);
    system->setStorePositions               (false);

    /*
    Optimizer* myOptimizer = new Optimizer(system, alpha, beta);
    myOptimizer->optimizeParameters();
    */

    system->runMetropolisSteps              (numberOfSteps);




/*
  Optimized parameters
  alpha : 0.9429569512
  beta  : 0.4603492076
*/

/*
Many
 Energy          : 3.0112
 Variance        : 2.9009e-09
 Acceptance rate : 0.999187




Two
 Energy          : 3.0141
 Variance        : 2.8288e-09
 Acceptance rate : 0.961014







*/



    return 0;
}





