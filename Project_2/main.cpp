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
    int numberOfSteps       = (int) 1e7;
    double omegaHO          = 1;           // Oscillator frequency.
    double omegaZ           = 1;
    double alpha            = 1;      // Variational parameter.
    double beta             = 0.3;          // Variational parameter.
    double gamma            = 2.82843;
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.1;          // Fraction steps used for equilibration.
    double C                = 1.0;
    double a                = 1;

    System* system = new System();
    system->setInitialState                 (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    if (1){
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
    system->setPrintProgress                (true);
    bool   optimizing =                      true;


    if (optimizing){
        Optimizer* myOptimizer = new Optimizer(system, alpha, beta);
        myOptimizer->optimizeParameters();

        system->getWaveFunction()->setAlpha(myOptimizer->getAlpha());
        system->getWaveFunction()->setBeta(myOptimizer->getBeta());
    }

    system->setStorePositions               (false);
    system->setPrintResults                 (true);
    system->runMetropolisSteps              (numberOfSteps);



 /* N2NoJ
 Alpha :      0.9507004936
 Beta  :      0.3
 omega :      1
 */

 /* N2HO
 Alpha :      1.003164596
 Beta  :      0.3
 omega :      1
 */




    /*
  Optimized parameters

  -- System info --
 Name : Many body quantum dot
 Number of particles  : 2
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 4
 Alpha :      1.00338
 Beta  :      0.3
 Omega :      1
 a     :      1

  ----- Reults -----
 Energy          : 2.9977
 Variance        : 1.6587e-08
 Acceptance rate : 0.999214






  -- System info --
 Name : Many body quantum dot
 Number of particles  : 6
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 4
 Alpha :      1.026740104
 Beta  :      0.4
 Omega :      1
 a     :      1

  ----- Reults -----
 Energy          : 20.19
 Variance        : 2.1611e-07
 Acceptance rate : 0.998167

*/



 /*
  -- System info --
 Name : Many body quantum dot
 Number of particles  : 12
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 4
 Alpha :      0.9
 Beta  :      0.5
 Omega :      1
 a     :      1

  ----- Reults -----
 Energy          : 65.878
 Variance        : 7.4929e-07
 Acceptance rate : 0.996098

*/


/*
  -- System info --
 Name : Many body quantum dot
 Number of particles  : 20
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 4
 Alpha :      0.9
 Beta  :      0.5
 Omega :      1
 a     :      1

  ----- Reults -----
 Energy          : 155.9
 Variance        : 1.7264e-06
 Acceptance rate : 0.992662
*/


/*
  -- System info --
 Name : Two body quantum dot
 Number of particles  : 2
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 5
 Alpha :      0.97
 Beta  :      0.43
 Omega :      1
 a     :      1

  ----- Reults -----
 Energy          : 3.0006
 Variance        : 2.0239e-09
 Acceptance rate : 0.96536
*/




/*
  -- System info --
 Name : Many body quantum dot
 Number of particles  : 20
 Number of dimensions : 2
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 4
 Alpha :      6.6
 Beta  :      0.2
 Omega :      0.1
 a     :      0.1

  ----- Reults -----
 Energy          : 30.51
 Variance        : 4.4969e-08
 Acceptance rate : 0.999656
 Mean distance   : 106.69
*/





    return 0;
}





