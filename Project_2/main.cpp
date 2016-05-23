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




    //    ofstream energyFile;
    //    energyFile.open("dataFiles/energyPlotfileN6.txt", ios::out);
    //    int na = 20;
    //    int nb = 20;

    //    for (int ai=0; ai<=na; ai++){
    //        for (int bi=0; bi<=nb; bi++){

    //            cout << "ai: "<< ai << "    bi: " << bi << "\r";
    //            fflush(stdout);

    int numberOfParticles   = 2;
    int numberOfDimensions  = 2;
    int numberOfSteps       = (int) 1e7;
    double omegaHO          = 1.;           // Oscillator frequency.
    double omegaZ           = 1.0;
    double alpha            = 1.00;//338;//0.94295;      // Variational parameter.
    double beta             = 0.3;      // Variational parameter.
    double gamma            = 2.82843;
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.1;          // Fraction steps used for equilibration.
    double C                = 1.0;
    double a                = 0;

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
    system->setStorePositions               (true);
    bool   optimizing =                      false;


    if (optimizing){
        Optimizer* myOptimizer = new Optimizer(system, alpha, beta);
        myOptimizer->optimizeParameters();

        system->getWaveFunction()->setAlpha(myOptimizer->getAlpha());
        system->getWaveFunction()->setBeta(myOptimizer->getBeta());
    }

    system->setPrintResults                 (true);
    system->runMetropolisSteps              (numberOfSteps);
    //energyFile << system->getSampler()->getEnergy() << "    ";
    //        }
    //        energyFile << endl;
    //    }
    //    energyFile.close();




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





    return 0;
}





