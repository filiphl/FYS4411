#include "manybodyquantumdotwavefunction.h"

ManyBodyQuantumDotWaveFunction::ManyBodyQuantumDotWaveFunction(System *system, double omega, double a, double alpha, double beta) :
    WaveFunction(system)
{
    m_omega = omega;
    m_a = a;
    m_alpha = alpha;
    m_beta = beta;
    m_npHalf =  m_system->getNumberOfParticles()/2;
    setupSlater();
}

double ManyBodyQuantumDotWaveFunction::computeRatio(std::vector<class Particle*> particles, int i, int j, double change)
{
    particles[i]->adjustNewPosition(change, j);
    double RSD = 0;
    double RC = 0;


    if (i<m_npHalf){
        for (int k=0; k<m_npHalf; k++){
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            RSD += m_slaterUpInverse(k,i)*computeSingleParticleWF(nx, ny, particles[i]->getNewPosition()[0], particles[i]->getNewPosition()[1]);
        }
    }
    else{
        for (int k=m_npHalf; k<m_npHalf*2; k++){
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            RSD += m_slaterUpInverse(k,i)*computeSingleParticleWF(nx, ny, particles[i]->getNewPosition()[0], particles[i]->getNewPosition()[1]);
        }
    }

    double exponentN = 0;
    double exponentO = 0;
    for(int k=0; k<m_npHalf*2-1; k++){
        for (int p=k+1; p<m_npHalf; p++){
            double r_kpN = 0;
            double r_kpO = 0;
            for (int d=0; d<m_system->getNumberOfDimensions(); d++){
                r_kpN += (particles[k]->getNewPosition()[d]-particles[p]->getNewPosition()[d])
                        *(particles[k]->getNewPosition()[d]-particles[p]->getNewPosition()[d]);
                r_kpO += (particles[k]->getOldPosition()[d]-particles[p]->getOldPosition()[d])
                        *(particles[k]->getOldPosition()[d]-particles[p]->getOldPosition()[d]);
            }
            r_kpN = sqrt(r_kpN);
            r_kpO = sqrt(r_kpO);
            exponentN += m_a*r_kpN/(1+m_beta*r_kpN);
            exponentO += m_a*r_kpO/(1+m_beta*r_kpO);
        }
    }
    RC = exp(exponentN - exponentO);
    return RSD*RC;
}

double ManyBodyQuantumDotWaveFunction::computeLaplacian(std::vector<class Particle*> particles) {

}

double ManyBodyQuantumDotWaveFunction::computeGradient(std::vector<Particle *> particles, int particle, int dimension){

}

double ManyBodyQuantumDotWaveFunction::evaluate(std::vector<class Particle*> particles)
{
    //double phii = computeSingleParticleWF(nx,ny, particles[i]->getOldPosition()[0], particles[i]->getOldPosition()[1] );
}


double ManyBodyQuantumDotWaveFunction::computeSingleParticleWF(int nx, int ny, double x, double y)
{
    return hermite(nx, x)*(sqrt(m_omega)*x) *
            hermite(ny, y)*(sqrt(m_omega)*y) *
            exp(-m_omega*(x*x + y*y)*0.5);
}




double ManyBodyQuantumDotWaveFunction::hermite(int energyLevel, double position)
{
    if (energyLevel == 0){
        return 1;
    }
    else if (energyLevel == 1) {
        return 2*position;
    }
    else if (energyLevel == 2) {
        return 4*position*position - 2;
    }
    else if (energyLevel == 3){
        return 8*position*position*position - 12*position;
    }
    else if (energyLevel == 4){
        return 16*position*position*position*position - 48*position*position + 12;
    }
    else {
        cout << "Too high energy level!"<<endl;
        exit(0);
    }
}


void ManyBodyQuantumDotWaveFunction::setupSlater()
{
    m_slaterUp   = zeros<mat>(m_npHalf, m_npHalf);
    m_slaterDown = zeros<mat>(m_npHalf, m_npHalf);

    mat m_quantumNumbers = zeros<mat>(10,2);
    m_quantumNumbers(0 ,0) = 0;  m_quantumNumbers(0, 1) = 0;
    m_quantumNumbers(1 ,0) = 1;  m_quantumNumbers(1, 1) = 0;
    m_quantumNumbers(2 ,0) = 0;  m_quantumNumbers(2, 1) = 1;
    m_quantumNumbers(3 ,0) = 2;  m_quantumNumbers(3, 1) = 0;
    m_quantumNumbers(4 ,0) = 1;  m_quantumNumbers(4, 1) = 1;
    m_quantumNumbers(5 ,0) = 0;  m_quantumNumbers(5, 1) = 2;
    m_quantumNumbers(6 ,0) = 3;  m_quantumNumbers(6, 1) = 0;
    m_quantumNumbers(7 ,0) = 2;  m_quantumNumbers(7, 1) = 1;
    m_quantumNumbers(8 ,0) = 1;  m_quantumNumbers(8, 1) = 2;
    m_quantumNumbers(9 ,0) = 0;  m_quantumNumbers(9, 1) = 3;

    for (int i=0; i<m_npHalf; i++){
        for (int j=0; j<m_npHalf; j++){
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            //cout << "nx=" << nx << "  ny=" << ny <<endl;
            m_slaterUp(i,j) = computeSingleParticleWF(nx,ny,m_system->getParticles()[i]->getOldPosition()[0], m_system->getParticles()[i]->getOldPosition()[1]);
            m_slaterDown(i,j) = computeSingleParticleWF(nx,ny,m_system->getParticles()[i+m_npHalf]->getOldPosition()[0], m_system->getParticles()[i+m_npHalf]->getOldPosition()[1]);
        }
    }




    m_slaterUpInverse   = m_slaterUp.i();
    m_slaterDownInverse = m_slaterDown.i();
    cout << "Initial slater up:\n"         << m_slaterUp << endl;
    cout << "Initial slater up inverse:\n" << m_slaterUpInverse << endl;
    cout << "Initial slater down:\n"         << m_slaterDown << endl;
    cout << "Initial slater down inverse:\n" << m_slaterDownInverse << endl;

}

void ManyBodyQuantumDotWaveFunction::updateSlater(int k)
{
    for (int j=0; j<m_npHalf; j++){
        int nx = m_quantumNumbers(j,0);
        int ny = m_quantumNumbers(j,1);
        m_slaterUp(k,j) = computeSingleParticleWF(nx, ny, m_system->getParticles()[k]->getNewPosition()[0], m_system->getParticles()[k]->getNewPosition()[1]);
    }
}



