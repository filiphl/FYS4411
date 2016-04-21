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

bool System::getAnalyticalLaplacian() const
{
    return m_AnalyticalLaplacian;
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

double System::getDerivativeStep() const
{
    return m_derivativeStep;
}

void System::setDerivativeStep(double derivativeStep)
{
    m_derivativeStep = derivativeStep;
}

void System::setParticles(const std::vector<Particle *> &particles)
{
    m_particles = particles;
}

bool System::metropolisStep() {

    int p = Random::nextInt(m_numberOfParticles);                   // Random particle
    int d = Random::nextInt(m_numberOfDimensions);                  // Random dimension

    if (m_importanceSampling){
        double qForceOld = qForce(p,d);
        prob = m_waveFunction->computeRatio(m_particles, p,d,dx);
        double qForceNew = qForce(p,d);
        double greensFunction = exp(0.5*(qForceOld+qForceNew)*((m_D*m_dt/2)*(qForceNew-qForceOld) - dx));           // Only term special for i,j = p,d

        prob *= greensFunction;
    }

    else{
        dx = m_stepLength*Random::nextGaussian(0,1/sqrt(2));
        prob = m_waveFunction->computeRatio(m_particles, p, d, dx);
    }


    double mynt = Random::nextDouble();              // Uniform [0,1]

    if (mynt < prob){ return true; }                 // Accept.

    else {                                           // Reject.
        m_particles[p]->adjustOldPosition(-dx, d);
        return false;
    }
}




void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < m_numberOfMetropolisSteps; i++) {

       /*if (i%100){     // Added by us.
            cout << "  " << setprecision(2) << 100*i/m_numberOfMetropolisSteps << "% complete"<< "\r";
            fflush(stdout);
        }
*/
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

    double force = 2*m_waveFunction->computeGradient(m_particles, i, j);

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
    sprintf(cmd, "rm %s", m_oldPositionFileName);
    system(cmd);
    m_oldPositionFile.open(m_oldPositionFileName, ios::out);
}


void System::closeEnergyFile()
{
    m_energyFile.close();
}


void System::closePositionFile()
{
    m_oldPositionFile.close();
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

void System::setAnalyticalLaplacian(bool value)
{
    m_AnalyticalLaplacian = value;
}




