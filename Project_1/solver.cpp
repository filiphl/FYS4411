#include "solver.h"
using namespace arma;

Solver::Solver(){}

void Solver::addparticle(){
    m_nParticles += 1;
}

double Solver::localenergy(mat r){
double kinetic = 0;
double potential = 0;
double hbar = 1;
double dstep = 0.1;
mat rplus = r;
mat rminus = r;
double psiplus;
double psiminus;
double psi = wavefunction(r);
for (int p = 0; p<m_nParticles; p++){
    double r_single = 0;
    for (int d=0; d<m_nDimensions; d++){
        rplus(p,d) = r(p,d) + dstep;
        rminus(p,d) = r(p,d) - dstep;

        psiplus = wavefunction(rplus);
        psiminus = wavefunction(rminus);


        kinetic += (psiplus -2*psi + psiminus);

        r_single += r(p,d)*r(p,d);

        rplus(p,d) = r(p,d);
        rminus(p,d) = r(p,d);
    }
    potential += r_single;
}
kinetic *= -0.5*(hbar*hbar)/(m_m*dstep*dstep*psi);
potential *= 0.5*m_m*m_w*m_w/psi;

return kinetic +potential;
}

double Solver::wavefunction(mat r){
    double psi = 1;
    for (int p=0; p<m_nParticles; p++){
        double r_single = 0;
        double g = 0;
        for (int d=0; d<m_nDimensions; d++){
            r_single += r(p,d) * r(p,d);
        }
        g = exp(-m_alpha * r_single );
        psi = psi*g;
    }
    return psi;
}


