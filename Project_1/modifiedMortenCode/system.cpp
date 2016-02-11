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
    double qForceOld = qForce(m_particles);

    double dx = m_stepLength * Random::nextGaussian(0,sqrt(m_dt)) + m_D*qForceOld*m_dt;  // sqrt(2*m_D*m_dt), but m_D=0.5.
    m_particles[p]->adjustPosition(dx , d);                         // Propose move
    double newWaveFunction = m_waveFunction->evaluate(m_particles);
    double qForceNew = qForce(m_particles);

    double greensFunction = 0.0;
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j< m_numberOfDimensions; j++){

            greensFunction += 0.5*(qForceOld+qForceNew)*
                    (m_D*m_dt*0.5*(-qForceOld+qForceNew)-m_particles[i]->getPosition()[j] + dx);
        }
    }
    greensFunction = exp(greensFunction);

    double prob = greensFunction * newWaveFunction*newWaveFunction / (oldWaveFunction*oldWaveFunction);
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


double System::qForce(std::vector<class Particle*> particles){
    double force = 0;
    for (int i=0; i<m_numberOfParticles; i++){
        for (int j=0; j<m_numberOfDimensions; j++){
            particles[i]->adjustPosition(-m_stepLength, j);
            double waveFunctionOld = m_waveFunction->evaluate(particles);
            particles[i]->adjustPosition(2*m_stepLength, j);
            double waveFunctionNew = m_waveFunction->evaluate(particles);
            particles[i]->adjustPosition(-m_stepLength, j);
            force += (waveFunctionNew - waveFunctionOld) /
                    (m_stepLength * m_waveFunction->evaluate(particles));
        }
    }
    return force;
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




