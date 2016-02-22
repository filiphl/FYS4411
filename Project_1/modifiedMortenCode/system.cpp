#include "system.h"


using namespace std;

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    int p = Random::nextInt(m_numberOfParticles);                   // Random particle
    int d = Random::nextInt(m_numberOfDimensions);                  // Random dimension
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);
    //double qForceOld = qForce(p,d);

    double dx = m_stepLength * Random::nextGaussian(0, sqrt(m_dt));// + m_D*qForceOld*m_dt;  // sqrt(2*m_D*m_dt), but m_D=0.5.
    m_particles[p]->adjustPosition(dx , d);                         // Propose move
    double newWaveFunction = m_waveFunction->evaluate(m_particles);
    //double qForceNew = qForce(p,d);

    //double greensFunction = -0.5*(qForceOld+qForceNew)*dx;           // Only term special for i,j = p,d
    /*
     for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j< m_numberOfDimensions; j++){

            m_particles[p]->adjustPosition(-dx , d);
            qForceOld = qForce(i, j);
            m_particles[p]->adjustPosition(dx , d);
            qForceNew = qForce(i, j);

            greensFunction += 0.5*(qForceOld+qForceNew) * m_D*m_dt*0.5*(qForceNew-qForceOld);
        }
    }
    */
    //greensFunction = exp(greensFunction);

    double prob = newWaveFunction*newWaveFunction / (oldWaveFunction*oldWaveFunction);//greensFunction *
    double mynt = Random::nextDouble();                             // Uniform [0,1]

    if (mynt < prob){   // Accept. Keep new position.
        return true;
    }

    else{               // Reject. Reset the position.
        m_particles[p]->adjustPosition(-dx, d);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < m_numberOfMetropolisSteps; i++) {

        if (i%100){     // Added by us.
            cout << "  " << setprecision(2) << 100*i/m_numberOfMetropolisSteps << "% complete"<< "\r";
            fflush(stdout);
        }

        bool acceptedStep = metropolisStep();
        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

        if (i > m_equilibrationFraction * m_numberOfMetropolisSteps) {
            m_sampler->sample(acceptedStep);
        }

    }
    m_sampler->computeAverages();
    m_sampler->computeAnalyticalEnergy();
    m_sampler->printOutputToTerminal();
}


double System::qForce(int i, int j){
    /*
    m_particles[i]->adjustPosition(-m_stepLength, j);               // -
    double waveFunctionOld = m_waveFunction->evaluate(m_particles);
    m_particles[i]->adjustPosition(2*m_stepLength, j);              // +
    double waveFunctionNew = m_waveFunction->evaluate(m_particles);
    m_particles[i]->adjustPosition(-m_stepLength, j);               // reset
    double force = (waveFunctionNew - waveFunctionOld) /
            (m_stepLength * m_waveFunction->evaluate(m_particles));
    */

    double force = -2*m_waveFunction->getParameters()[0]*m_particles[i]->getPosition()[j];
    return force;
}

void System::openFile()
{
    m_outfile.open(m_filename, ios::out);
}

void System::closeFile()
{
    m_outfile.close();
}




void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}




