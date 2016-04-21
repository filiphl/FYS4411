#include "manybodyquantumdotwavefunction.h"

ManyBodyQuantumDotWaveFunction::ManyBodyQuantumDotWaveFunction(System *system) :
    WaveFunction(system)
{
    setupSlater();
}

double ManyBodyQuantumDotWaveFunction::computeRatio(std::vector<class Particle*> particles, int i, int j, double change)
{

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
    cout << "in SingleParticleWf"<<endl;
    return hermite(nx, x)*(sqrt(m_omega)*x) *
           hermite(ny, y)*(sqrt(m_omega)*y) *
           exp(-m_omega*(x*x + y*y)*0.5);
}




double ManyBodyQuantumDotWaveFunction::hermite(int energyLevel, double position)
{
   if (energyLevel == 0){
       cout << "in hermite"<<endl;
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
    int npHalf = m_system->getNumberOfParticles()/2;

    m_slaterUp   = zeros<mat>(npHalf, npHalf);
    cout <<"slaterup: "<< m_slaterUp<<endl;
    m_slaterDown = zeros(npHalf, npHalf);

    mat quantumNumbers = zeros<mat>(10,2);
    quantumNumbers(0 ,0) = 0;  quantumNumbers(0, 1) = 0;
    quantumNumbers(1 ,0) = 1;  quantumNumbers(1, 1) = 0;
    quantumNumbers(2 ,0) = 0;  quantumNumbers(2, 1) = 1;
    quantumNumbers(3 ,0) = 2;  quantumNumbers(3, 1) = 0;
    quantumNumbers(4 ,0) = 1;  quantumNumbers(4, 1) = 1;
    quantumNumbers(5 ,0) = 0;  quantumNumbers(5, 1) = 2;
    quantumNumbers(6 ,0) = 3;  quantumNumbers(6, 1) = 0;
    quantumNumbers(7 ,0) = 2;  quantumNumbers(7, 1) = 1;
    quantumNumbers(8 ,0) = 1;  quantumNumbers(8, 1) = 2;
    quantumNumbers(9 ,0) = 0;  quantumNumbers(9, 1) = 3;


    cout << "npHalf = "<<npHalf<< endl;
    for (int i=0; i<npHalf; i++){
        cout << "i = "<< i <<endl;
        for (int j=0; j<npHalf; j++){
            cout << "j = "<< j <<endl;
            cout << m_system->getParticles()[i]->getOldPosition()[0]<<endl;
            m_slaterUp(i,j) = computeSingleParticleWF(i,j,m_system->getParticles()[i]->getOldPosition()[0], m_system->getParticles()[i]->getOldPosition()[1]);
            cout << "TEST3"<<endl;
        }
    }

    cout << m_slaterUp << endl;


}



