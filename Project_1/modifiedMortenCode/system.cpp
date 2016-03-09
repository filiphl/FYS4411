#include "system.h"


using namespace std;

void System::setImportanceSampling(bool value)
{
    m_importanceSampling = value;
}

void System::setStoreLocalEnergy(bool value)
{
    m_storeLocalEnergy = value;
    openEnergyFile();
}

bool System::getAnalyticalDoublederivative() const
{
    return m_analyticalDoublederivative;
}

bool System::getImportanceSampling() const
{
    return m_importanceSampling;
}

bool System::getStoreLocalEnergy() const
{
    return m_storeLocalEnergy;
}

bool System::OptimizingParameters() const
{
    return m_optimizingParameters;
}

void System::OptimizingParameters(bool optimizingParameters)
{
    m_optimizingParameters = optimizingParameters;
}

bool System::getStorePositions() const
{
    return m_storePositions;
}

void System::setStorePositions(bool storePositions)
{
    m_storePositions = storePositions;
    openPositionFile();
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    int p = Random::nextInt(m_numberOfParticles);                   // Random particle
    int d = Random::nextInt(m_numberOfDimensions);                  // Random dimension

    if (m_importanceSampling){
        double oldWaveFunction = m_waveFunction->evaluate(m_particles);
        double qForceOld = qForce(p,d);
        dx = m_stepLength*(Random::nextGaussian(0,sqrt(m_dt)) + m_D*qForceOld*m_dt);  // sqrt(2*m_D*m_dt), but m_D=0.5.
        m_particles[p]->adjustPosition(dx , d);                                        // Propose move
        double newWaveFunction = m_waveFunction->evaluate(m_particles);
        double qForceNew = qForce(p,d);
        double greensFunction = exp(0.5*(qForceOld+qForceNew)*((m_D*m_dt/2)*(qForceNew-qForceOld) - dx));           // Only term special for i,j = p,d

        prob = newWaveFunction*newWaveFunction / (oldWaveFunction*oldWaveFunction) * greensFunction;
    }

    else{
        double oldWaveFunction = m_waveFunction->evaluate(m_particles);
        dx = m_stepLength * (Random::nextDouble()*2-1);
        m_particles[p]->adjustPosition(dx , d);                          // Propose move
        //cout << "position: "<< m_particles[p]->getPosition()[d]<< endl;
        double newWaveFunction = m_waveFunction->evaluate(m_particles);

        prob = newWaveFunction*newWaveFunction / (oldWaveFunction*oldWaveFunction);
    }


    double mynt = Random::nextDouble();         // Uniform [0,1]

    if (mynt < prob){ return true; }            // Accept. Keep new position.

    else{                                       // Reject. Reset the position.
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
    if (m_storeLocalEnergy){ closeEnergyFile(); }
    if (m_storePositions){ closePositionFile(); }
}


double System::qForce(int i, int j){
/*
    m_particles[i]->adjustPosition(-m_derivativeStep, j);               // -
    double waveFunctionOld = m_waveFunction->evaluate(m_particles);
    m_particles[i]->adjustPosition(2*m_derivativeStep, j);              // +
    double waveFunctionNew = m_waveFunction->evaluate(m_particles);
    m_particles[i]->adjustPosition(-m_derivativeStep, j);               // reset
    double force = (waveFunctionNew - waveFunctionOld) /
            (m_derivativeStep * m_waveFunction->evaluate(m_particles)); // Factors of 2 cancel.
*/
    double force = -2*m_waveFunction->getParameters()[0]*m_particles[i]->getPosition()[j];
    return force;
}


void System::openEnergyFile()
{
    char cmd[50];
    sprintf(cmd, "rm %s", m_energyFileName);
    system(cmd);
    m_energyFile.open(m_energyFileName, ios::out);
}


void System::openPositionFile()
{
    char cmd[50];
    sprintf(cmd, "rm %s", m_positionFileName);
    system(cmd);
    m_positionFile.open(m_positionFileName, ios::out);
}


void System::closeEnergyFile()
{
    m_energyFile.close();
}


void System::closePositionFile()
{
    m_positionFile.close();
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

void System::setAnalyticalDoubleDerivative(bool value)
{
    m_analyticalDoublederivative = value;
}




