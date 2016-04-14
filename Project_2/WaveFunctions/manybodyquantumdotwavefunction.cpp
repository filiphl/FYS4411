#include "manybodyquantumdotwavefunction.h"
#include "particle.h"

ManyBodyQuantumDotWaveFunction::ManyBodyQuantumDotWaveFunction(System *system):
    WaveFunction(system)
{

}

double ManyBodyQuantumDotWaveFunction::evaluate(std::vector<Particle *> particles)
{

}

double ManyBodyQuantumDotWaveFunction::computeLaplacian(std::vector<Particle *> particles)
{

}

double ManyBodyQuantumDotWaveFunction::computeGradient(std::vector<Particle *> particles, int particle, int dimension)
{

}


/*
double ManyBodyQuantumDotWaveFunction::evaluate(std::vector<Particle *> particles)
{

    int energyLevel = nx+ny;
    double phii = computeSingleParticleWF(nx,ny, particles[i]->getPosition()[0], particles[i]->getPosition()[1] );
    return 0.0;

}


double ManyBodyQuantumDotWaveFunction::computeSingleParticleWF(int nx, int ny, double x, double y)
{

    return A *
           hermite(nx, x)*(sqrt(m_omega)*x) *
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
       return 0;
   }
}

*/

