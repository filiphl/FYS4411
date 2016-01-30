#include "solver.h"
using namespace arma;

Solver::Solver(){}

void Solver::addparticle(){
    m_nParticles += 1;
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
    potential *= 0.5*m_m*m_w*m_w;
    cout <<"Kinetic: " <<kinetic << "    Potetial: " << potential << endl;

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


mat Solver::metropolis_step(mat r){
    mat oldr = r;
    mat newr = r;
    double oldE;
    double newE;
    for (int p=0; p<m_nParticles; p++){
        for (int d=0; d<m_nDimensions; d++){
            newr(p,d) += Random::nextGaussian(0,1); //Random number ND
        }
        //cout << oldr << endl;
        oldE = localenergy(oldr);
        newE = localenergy(newr);
        if (newE<oldE){
            oldr = newr;}
        else{
            double mynt = Random::nextDouble();
            //cout << mynt <<endl;
            if (mynt < wavefunction(newr)/wavefunction(oldr)){
                oldr = newr;}
            else { newr = oldr;}
        }
    }
    return newr;


}



















