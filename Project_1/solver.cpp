#include "solver.h"
using namespace arma;

Solver::Solver(){
    m_dx = 1;
}

void Solver::addparticle(){

}


double Solver::Analytical(mat r){
    double energy = 0;
    double h2 = 1;
    for (int p=0; p<m_nParticles; p++){
        double r_single = 0;
        for (int d=0; d<m_nDimensions; d++){
            r_single += r(p,d)*r(p,d);
        }
        energy += h2*m_alpha*(-2*m_alpha*r_single + m_nDimensions) + 0.5*m_m*m_w*m_w*r_single;
    }
    return energy;
}


double Solver::localenergy(mat r){
    double kinetic = 0;
    double potential = 0;
    double hbar = 1;
    double dstep = 0.0001;
    mat rplus = r;
    mat rminus = r;
    double psiplus;
    double psiminus;
    double psi = wavefunction(r);


    for (int p = 0; p<m_nParticles; p++){
        double r_single = 0;
        for (int d=0; d<m_nDimensions; d++){
            rplus(p,d) += dstep;
            rminus(p,d) -= dstep;

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
    potential *= 0.5*m_m*m_w*m_w;


    return kinetic +potential;
}

double Solver::wavefunction(mat r){
    double psi = 1;
    for (int p=0; p<m_nParticles; p++){
        double r_single = 0;
        double g = 0;
        for (int d=0; d<m_nDimensions; d++){
            r_single += r(p,d)*r(p,d);
        }
        g = exp(-m_alpha * r_single );
        psi *= g;
    }
    return psi;
}


mat Solver::metropolis_step(mat r){
    mat oldr = r;
    mat newr = r;

    double oldwavefunction = wavefunction(oldr);
    for (int p=0; p<m_nParticles; p++){
        for (int d=0; d<m_nDimensions; d++){
            newr(p,d) += m_dx*Random::nextGaussian(0,0.2); //Random number ND
        }

        double newwavefunction = wavefunction(newr);
        double prob = newwavefunction*newwavefunction/(oldwavefunction*oldwavefunction);
        double mynt = Random::nextDouble();             // Uniform [0,1]

        if (mynt < prob){
            for (int d=0; d<m_nDimensions; d++){
                oldr(p,d) = newr(p,d);
            }
            oldwavefunction = newwavefunction;
            m_accepted++;
        }
        else{
            for (int d=0; d<m_nDimensions; d++){
                newr(p,d) = oldr(p,d);
            }
        }
    }

    return newr;
}

