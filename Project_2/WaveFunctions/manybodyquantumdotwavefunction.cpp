#include "manybodyquantumdotwavefunction.h"

ManyBodyQuantumDotWaveFunction::ManyBodyQuantumDotWaveFunction(System *system):
    WaveFunction(system)
{

}

double ManyBodyQuantumDotWaveFunction::evaluate(std::vector<Particle *> particles, int nx, int ny)
{
    int energyLevel = nx+ny;
    double phii = computeSingleParticleWF(nx,ny, particles[i]->getPosition()[0], particles[i]->getPosition()[1] )
}


double ManyBodyQuantumDotWaveFunction::computeSingleParticleWF(int nx, int ny, double x, double y)
{
    return A * hermite(nx, x)*(sqrt(m_omega)*x) * hermite(ny, y)*(sqrt(m_omega)*y) * exp(-m_omega*(x*x + y*y)*0.5);
}




double ManyBodyQuantumDotWaveFunction::hermite(int energyLevel, double position)
{
   if (energyLevel == 0){
        return H0;
   }
   else if (energyLevel == 1) {
        return H1;
   }
   else if (energyLevel == 2) {
        return H2;
   }
   else if (energyLevel == 3){
        return H3;
   }
   else {
       cout << "Too high energy level!"<<endl;
       exit(0);
   }
}



